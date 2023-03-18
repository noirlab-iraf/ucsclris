# MAPMASK: Map RA, DEC into X,Y on slitmask
#
# This routine contains a conversion of the fortran code lifted from AUTOSLIT.
# It has been slightly modified to reference all positions relative to
# the field center, which is identified as the center of the mask.
#
# For simplicity, we read the input file once; write a temporary _sorted_ file
# and then read this again to work on adjusting the slit lengths.

# Rev:  3-jul-97 --- put in sep_slit to describe a setable spacing between slits
# Rev: 10-jul-97 --- add TV pixels for Guide Star.
# Rev: 12-aug-97 --- add slit-ccd mapping; added update to PA for precession
# 26mar98: Allow PCODE=0 objects as reference -- carry through, but ignore
# 17apr98: Replaced map_sort with more robust version
# 12jan07: Add ADC switch

include	<math.h>
include	<mach.h>
include	<error.h>
include <math/gsurfit.h>

define	LAT	19.8265		# Keck latitude
define	M_BEND	1.94		# Mask bend in degrees
define	M_ANGLE	8.06		# Mask angle (tilt) in degrees
define	MASK_X0	107.95		# Half of Mask X size (215.9) (mm)
define	MASK_Y0	167.8		# Half of Mask Y size (335.6) (mm)
define	MM_ARCS	0.7253		# mm per arcsec
define	ADC_ADJ	0.99857d0	# plate scale adjustment for ADC
define	MM_PIX	0.15767		# mm per pixel
define	FP_XCEN	305.		# Mask x-center (mm) in FP coord
define	FP_YCEN	0.		# Mask y-center (mm) in FP coord
define	CCDXCEN	1024.5		# CCD x-center
define	CCDYCEN	1024.5		# CCD y-center
define	XPIXOFF	0.		# x-offset(pix) (not in AUTOSLIT)
define	YPIXOFF	80. 		# y-offset(pix)
define	MIN_MY	3.5		# Minimum mask y-value
define	MAX_MY	332.		# Maximum mask y-value
define	YSCLCOR	1.0175		# yscale_cor
define	SCL_FIX	1.0029		# Scaling correction
define	ACORR1	-2.95		# Astrom. correction term 1
define	ACORR2	7.734e-3	# Astrom. correction term 2
define	ACORR3	-3.212e-6	# Astrom. correction term 3
define	ACORR4	-3.274e-10	# Astrom. correction term 4

define	TV_XOFF		231.8 	# edge of TV CCD (mm, FPCS)
define	TV_YOFF		31.5	# edge of TV CCD (mm, FPCS)
define	TV_MMPIX	0.1835	# pixel scale: mm (FP) / pix (TV)
#define	TV_XOFF		233.5	# edge of TV CCD (mm, FPCS)
#define	TV_YOFF		30.5	# edge of TV CCD (mm, FPCS)

define		SZ_INLN	128	# Maximum length of input line
define		SZ_ID	32	# Maximum length of ID string
define		CODE_RF	0	# Priority code for reference object	XXX
define		CODE_GS	-1	# Priority code for guide star		XXX
define		CODE_AS	-2	# Priority code for alignment star	XXX

procedure	t_mapmask()

char	objfile[SZ_FNAME]			# ID, prior, mag, RA, Dec
char	output[SZ_FNAME]			# output file name
real	pa_mask
double	ha0					# Hour angle (field center)
double	pres, wave0, temp
double	epoch					# Epoch for PA precession
real	sep_slit				# minimum separation of slits
real	def_ybot, def_ytop, box_rad		# default (min) sizes
bool	adc					# ADC in?
bool	verbose					# verbose output?
real	prescan					# columns of CCD prescan
pointer	fda, fdb				# in/out file descriptors

int	pcode				# formerly "priority"
double	equinox				# copied value for coord. equinox
real	magn, pa, ybot, ytop
real	delpa				# for precessing PA

# char	sortfile[SZ_FNAME]		# temporary sort file
# pointer	fdx

char	tchar, id[SZ_ID]
int	npts, ndx, i, j
int	nslits				# number of actual slits
real	y1, y2, ytot, ylow, yupp, yb, yt	# used for verbose output
double	guide_ra, guide_dec		# saved coords of guide star
pointer	bufid, bufndx, bufmag
pointer	bufx, bufy, bufy1, bufy2, bufpc, bufpa, bufyb, bufyt
pointer	bufxccd, bufyccd		# ccd vectors

double	cosa, sina
double	ra0, dec0, ra, dec		# field center, object
double	ha				# HA of object
double	rlat				# latitude in radians
double	xaxis				# tot. offset (arcsec) to tel.
double	ra_ref, dec_ref			# refracted coords
double	delra_tel, deldec_tel		# apparent offsets (rad) to tel.
double	ra_tel, dec_tel			# apparent RA, dec (rad) of tel.

double	xfp, yfp			# "ideal" coords in focal plane system
double	x0, y0				# calculated x,y of field center
double	delx, dely			# calculated delta x,y from field center

real	x, y

int	guide_star			# is there a guide star?

int	xtv, ytv			# x,y of guide star in TV pixels
real	xgs, ygs			# x,y of guide star in FPCS (mm)

bool	clgetb(), strne()
real	clgetr()
int	map_sort()
int	fscan(), nscan()
pointer	open()

begin
	call clgstr ("objfile", objfile, SZ_FNAME)
	call clgstr ("output", output, SZ_FNAME)
	ha0 = clgetr ("ha0")
	adc = clgetb ("ADC")
	temp = clgetr ("temp")
	pres = clgetr ("pressure")
	wave0 = clgetr ("lambda_cen") * 1.e-4		# in microns
	sep_slit = clgetr ("sep_slit") * MM_ARCS	# sep. bet. slits (mm)
# (don't bother to correct for ADC scale, too small ...)
	def_ybot = clgetr ("lower_min")			# min. length in arcsec
	def_ytop = clgetr ("upper_min")			# min. length in arcsec
	box_rad = 0.5 * clgetr ("box_sz")		# box 1/2-length in arcs
	prescan = clgetr ("prescan")
	verbose = clgetb ("verbose")
	epoch = clgetr ("epoch")

        fda = open (objfile, READ_ONLY, TEXT_FILE)

# Count the entries
	ndx = 0
	while (fscan (fda) != EOF) {
		call gargwrd (tchar, 1)
		if (tchar == '#' || nscan() == 0) {
			next
		}
		call reset_scan()
		ndx = ndx + 1
	}
	npts = ndx
	call seek (fda, BOF)

# Allocate arrays
	call malloc (bufx, npts, TY_REAL)
	call malloc (bufy, npts, TY_REAL)
	call malloc (bufy1, npts, TY_REAL)
	call malloc (bufy2, npts, TY_REAL)
	call malloc (bufpc, npts, TY_INT)
	call malloc (bufpa, npts, TY_REAL)
	call malloc (bufyb, npts, TY_REAL)
	call malloc (bufyt, npts, TY_REAL)
	call malloc (bufmag, npts, TY_REAL)
	call malloc (bufxccd, npts, TY_REAL)
	call malloc (bufyccd, npts, TY_REAL)
	call malloc (bufid, npts*SZ_ID, TY_CHAR)
	call amovkc (EOS, Memc[bufid], npts*SZ_ID)
	call malloc (bufndx, npts, TY_INT)

# Now get Field Center and work out its lris coords
	while (fscan(fda) != EOF) {
		call gargwrd (tchar, 1)
		if (tchar == '#' || nscan() == 0) {
			next
		}
		call reset_scan()
		call gargwrd (id, SZ_ID)
		call gargi (pcode)
		call gargd (equinox)
		call gargd (ra0)
		call gargd (dec0)
		call gargr (pa_mask)
		if (nscan() < 6 || strne (id, "CENTER"))
		    call fatal (0, "CENTER not first line or poorly formatted")

# Convert to radians
		ha0  = DEGTORAD (15. * ha0)
		ra0  = DEGTORAD (15. * ra0)
		dec0 = DEGTORAD (dec0)
		break
	}

	rlat = DEGTORAD(LAT)

# Work out apparent tel. center based on refracted coords of field center
	call rad_refract (ra0, dec0, ha0, rlat, pres, temp, wave0,
								ra_ref, dec_ref)
	cosa = cos (DEGTORAD(-pa_mask))
	sina = sin (DEGTORAD(-pa_mask))

# NEW: adc adjustment
	if (adc) {
		xaxis = FP_XCEN / (MM_ARCS * ADC_ADJ)
	} else {
		xaxis = FP_XCEN / MM_ARCS
	}
	delra_tel = DEGTORAD(xaxis/3600.) * cosa / cos (dec_ref)
	deldec_tel = DEGTORAD(xaxis/3600.) * sina

	ra_tel = ra_ref - delra_tel
	dec_tel = dec_ref - deldec_tel

# Here's the calculated field center
	call fp_coord (ra_ref, dec_ref, ra_tel, dec_tel, cosa, sina, xfp, yfp, adc)

	call lris_coord (xfp, yfp, 0., 0., x0, y0)

#
# Loop through the rest of the file, keeping (x,y,pr,pa,y1,y2) and input lines
	ndx = 0
	guide_star = NO
	while (fscan(fda) != EOF) {
		call gargwrd (tchar, 1)
		if (tchar == '#' || nscan() == 0) {
			next
		}
		call reset_scan()
		call gargwrd (id, SZ_ID)
		call gargi (pcode)
		call gargr (magn)
		call gargd (ra)
		call gargd (dec)
		call gargr (pa)
		call gargr (ybot)
		call gargr (ytop)
		if (nscan() < 5) {
		    call eprintf ("Poorly-formatted data line -- skipped\n")
		    next
		}

# If guide (pickoff mirror) star, save RA and Dec
		if (pcode == CODE_GS) {
			if (guide_star == YES) {
				call eprintf ("Additional Guide Star ignored\n")
				next
			}
			guide_star = YES
			guide_ra  = ra
			guide_dec = dec
		}

# Additional slit info to be carried along: individual PA ...
		if (nscan() < 6 || pa == INDEF)
			pa = pa_mask
		pa = mod ((pa - pa_mask + 360.), 180.)
		if (pa > 90.)
			pa = pa - 180.
# ... and lengths
		if (nscan() < 7 || ybot == INDEF)
			ybot = def_ybot
		if (nscan() < 8 || ytop == INDEF)
			ytop = def_ytop
# Check codes and force the right sizes:
		if (pcode == CODE_GS) {
			ybot = 0.5 / MM_ARCS		# give it 1.0mm
			ytop = 0.5 / MM_ARCS		# just a marker, no adj
		} else if (pcode == CODE_AS) {
			ybot = box_rad
			ytop = box_rad
		} else if (pcode == CODE_RF) {
			ybot = 0.			# for safety
			ytop = 0.
		}

# Ready to calculate the mask coords. First convert RA,Dec to radians
		ra  = DEGTORAD (15. * ra)
		dec = DEGTORAD (dec)

# refract the coordinates; work out focal plane, lris coords:
		ha = ha0 + (ra0 - ra)
		call rad_refract (ra, dec, ha, rlat, pres, temp, wave0,
								ra_ref, dec_ref)

		call fp_coord (ra_ref, dec_ref, ra_tel, dec_tel, cosa, sina,
								xfp, yfp, adc)

		call lris_coord (xfp, yfp, x0, y0, delx, dely)

		x = delx + MASK_X0			# not used; TMP?
		y = dely + MASK_Y0			# not used; TMP?

# Check to make sure it's on mask
		if ((y < MIN_MY || y > MAX_MY) && pcode != CODE_GS) {
			call eprintf ("Object %s off mask -- skipped\n")
				call pargstr (id, SZ_ID)
			next
		}

		Memi[bufndx+ndx] = ndx
		Memi[bufpc+ndx] = pcode
		Memr[bufx+ndx] = delx
		Memr[bufy+ndx] = dely
		Memr[bufpa+ndx] = pa
# NEW! adjust for ADC scale change (overkill)
		if (adc) {		
			Memr[bufyb+ndx] = ybot * (MM_ARCS * ADC_ADJ)
			Memr[bufyt+ndx] = ytop * (MM_ARCS * ADC_ADJ)
		} else {
			Memr[bufyb+ndx] = ybot * MM_ARCS
			Memr[bufyt+ndx] = ytop * MM_ARCS
		}
		call strcpy (id, Memc[bufid+ndx*SZ_ID], SZ_ID)
		Memr[bufmag+ndx] = magn
		ndx = ndx + 1
	}
	npts = ndx
	call close (fda)

## OK, we have the (x,y) values. Rewrite input sorted in y..................

	nslits = map_sort (Memr[bufy], npts, Memr[bufx], Memr[bufyb],
			Memr[bufyt], Memi[bufndx], Memi[bufpc], guide_star)

# Guide star (if present) now last in list; if present, save x,y FPCS
	if (Memi[bufpc+npts-1] == CODE_GS) {
		xgs = FP_XCEN - Memr[bufx+npts-1] * cos (DEGTORAD (M_ANGLE))
		ygs = FP_YCEN - Memr[bufy+npts-1]
	}


## Now ready to work out the slit lengths.............................
# Recall, these are based on delta values so we can treat the bar sensibly
# Only objects which fall on mask included in list, and the guide star -- if
# present -- is at end.

# First, assign nominal limits to slits:

	call asubr (Memr[bufy], Memr[bufyt], Memr[bufy2], npts)
	call aaddr (Memr[bufy], Memr[bufyb], Memr[bufy1], npts)

# ... then make the edges butt
	call slit_adjust (Memr[bufy], Memr[bufyb], Memr[bufyt], Memr[bufy1],
				Memr[bufy2], Memi[bufpc], nslits, sep_slit)


# If verbose, print out differences between request and actual lengths:
	if (verbose) {
		call printf ("\nLen('')  Len(mm)   bot.,top    (requested)     del1,del2   Code  Ident\n")

		do i = 0, nslits-1 {
			ndx = Memi[bufndx+i]
			y = Memr[bufy+i]
			y1 = Memr[bufy1+i]
			y2 = Memr[bufy2+i]
			yb = Memr[bufyb+i]
			yt = Memr[bufyt+i]
			ytot = y1 - y2
			ylow = y1 - y  + 4.9e-4		# last term for roundoff
			yupp = y  - y2 + 4.9e-4

			call printf (
	"%6.2f %8.2f %6.2f,%5.2f  (%5.2f,%5.2f)  %5.2f,%5.2f %1s %4d  %s\n")
# (don't bother to adjust for ADC, as too small for this text ...)
				call pargr (ytot / MM_ARCS)
				call pargr (ytot)
				call pargr (ylow)
				call pargr (yupp)
				call pargr (yb)
				call pargr (yt)
				call pargr (ylow - yb)
				call pargr (yupp - yt)
# ... work out status code
			    if (yb > ylow || yt > yupp) {
				if ((ylow + yupp) >= (yt + yb))
					call pargc ("^")	# Size OK
				else
					call pargc ("*")	# Too short
			    } else {
				call pargc (" ")		# All OK
			    }
			    call pargi (Memi[bufpc+i])
			    call pargstr (Memc[bufid+ndx*SZ_ID], SZ_ID)
		}
	}

# ... Now update the coords to the center of the mask
	call aaddkr (Memr[bufx],  MASK_X0, Memr[bufx],  npts)
	call aaddkr (Memr[bufy],  MASK_Y0, Memr[bufy],  npts)
	call aaddkr (Memr[bufy1], MASK_Y0, Memr[bufy1], npts)
	call aaddkr (Memr[bufy2], MASK_Y0, Memr[bufy2], npts)

# Find appropriate slit limits for CODE_RF objects:
	do i = nslits, npts-1 {
		if (Memi[bufpc+i] == CODE_RF) {
			y = Memr[bufy+i]
			do j = 0, nslits-1 {
				if (Memr[bufy1+j] > y && Memr[bufy2+j] < y) {
					Memr[bufy1+i] = Memr[bufy1+j]
					Memr[bufy2+i] = Memr[bufy2+j]
					break
				}
			}
		}
	}


# Now work out predicted CCD coordinates
	call ccd_map (Memr[bufx], Memr[bufy], Memr[bufxccd], Memr[bufyccd], npts)

# One final kludge: precession will cause a change in PA -- until precession
# is worked into the whole package, calculate the updated PA for now.
# Precession approximation from Lang, using m,n for 1975.0
	delpa = RADTODEG (atan ((epoch - equinox) * 9.7157e-5 * sin (ra0) / cos (dec0)))

### Finally ready to assemble output file
# Open and write header info
        fdb = open (output, NEW_FILE, TEXT_FILE)

	call fprintf (fdb,
	"###########\n##\n##     MAPMASK:   Input = %s\n##\n###########\n##\n")
		call pargstr (objfile, SZ_FNAME)

	call fprintf (fdb,
	"##     Lambda (A)     T (C)      P (mm Hg)      HA (hr)    Prescan\n")
	call fprintf (fdb,
	"##     %7.1f        %4.1f      %7.1f          %4.1f   %8.0f\n")  
		call pargd (wave0 * 1.e4)
		call pargd (temp)
		call pargd (pres)
		call pargd (RADTODEG(ha0)/15.)
		call pargr (prescan)

# Make sure that dec is properly formatted for Keck target list
	if (dec0 < 0.) {
		call fprintf (fdb,
			"##\n##  Field_Center.... %011.2h  -%010.1h  %7.2f\n")
	} else {
		call fprintf (fdb,
			"##\n##  Field_Center.... %011.2h   %010.1h  %7.2f\n")
	}
		call pargd (RADTODEG(ra0) / 15.)
		call pargd (abs (RADTODEG(dec0)))
		call pargd (equinox)

	if (guide_star == YES) {
	    if (guide_dec < 0.) {
		call fprintf (fdb,
			"##\n##  Guide_Star...... %011.2h  -%010.1h  %7.2f\n")
	    } else {
		call fprintf (fdb,
			"##\n##  Guide_Star...... %011.2h   %010.1h  %7.2f\n")
	    }
		call pargd (guide_ra)
		call pargd (abs (guide_dec))
		call pargd (equinox)

		xtv = (TV_YOFF - ygs) / TV_MMPIX + 0.5
		ytv = (xgs - TV_XOFF) / TV_MMPIX + 0.5
		call fprintf (fdb,
		    "##\n##  Guide Star, Slit-Viewing TV: %dx %dy \n")
			call pargi (xtv)
			call pargi (ytv)
	} else {
		call fprintf (fdb, "##\n##  Guide_Star...... (none)\n")
	}

	call fprintf (fdb,
	"##\n##     Mask Position Angle: %8.3f deg in %7.2f   (%6.1f in %6.1f)\n##\n")
		call pargr (pa_mask + delpa)
		call pargd (epoch)
		call pargr (pa_mask)
		call pargd (equinox)

	call fprintf (fdb, "##\n## Xobj   Yobj     Ymin    Ymax   Rel_PA  CCDx   CCDy  Code  Mag  ID\n##\n")

# Loop through objects/slit and print:
	do i = 0, npts-1 {
		ndx = Memi[bufndx+i]
		call fprintf (fdb,
		    "%7.3f %7.3f %8.3f %7.3f %6.2f %7.1f %6.1f %4d %5.2f %s\n")
			call pargr (Memr[bufx+i])
			call pargr (Memr[bufy+i])
			call pargr (Memr[bufy1+i])
			call pargr (Memr[bufy2+i])
			call pargr (Memr[bufpa+ndx])
			call pargr (Memr[bufxccd+i])
			call pargr (Memr[bufyccd+i])
			call pargi (Memi[bufpc+i])
			call pargr (Memr[bufmag+ndx])
			call pargstr (Memc[bufid+ndx*SZ_ID], SZ_ID)
	}

	call close (fdb)

	call mfree (bufndx, TY_INT)
	call mfree (bufid, TY_CHAR)
	call mfree (bufpc, TY_INT)
	call mfree (bufyccd, TY_REAL)
	call mfree (bufxccd, TY_REAL)
	call mfree (bufmag, TY_REAL)
	call mfree (bufyt, TY_REAL)
	call mfree (bufyb, TY_REAL)
	call mfree (bufpa, TY_REAL)
	call mfree (bufpc, TY_REAL)
	call mfree (bufy2, TY_REAL)
	call mfree (bufy1, TY_REAL)
	call mfree (bufy, TY_REAL)
	call mfree (bufx, TY_REAL)

end


# CCD_MAP: map x,y (slitmask) into CCD coords. Assumes fit from geomap entered
# in def_map()

# To produce the coeffs, run geomap on a file that has X,Y,x,y (SM,CCD).
# The two columns of coeffs are for the x-fit and y-fit respectively.
# Note that the linear terms (surface1) must be added to the first 3
# coeffs of the distortion terms (surface2)

procedure	ccd_map (xmm, ymm, xccd, yccd, npts)

real	xmm[npts], ymm[npts]			# input slitmask coords
real	xccd[npts], yccd[npts]			# output ccd coords
int	npts					# number of slits

pointer	sfx, sfy

begin
	call def_ccd_map (sfx, sfy)
	call gsvector (sfx, xmm, ymm, xccd, npts)
	call gsvector (sfy, xmm, ymm, yccd, npts)
	call gsfree (sfx)
	call gsfree (sfy)

end

# DEF_CCD_MAP: define the slit-ccd mapping

procedure	def_ccd_map (sfx, sfy)

# These coefficients were taken from a fit to the M71_D1 map.
# Additional offsets correct to LRIS stow position (zeropt for flexure mapping)

pointer	sfx, sfy			# pointers to surface fits in x,y

int	ncoeff
pointer	xcoeff, ycoeff			# coeff's in x,y

begin
	ncoeff = 20
	call malloc (xcoeff, ncoeff, TY_REAL)
	call malloc (ycoeff, ncoeff, TY_REAL)

	Memr[xcoeff  ]  = 2.
	Memr[xcoeff+1]  = 4.
	Memr[xcoeff+2]  = 3.
	Memr[xcoeff+3]  = 1.
	Memr[xcoeff+4]  = 1.
	Memr[xcoeff+5]  = 216.
	Memr[xcoeff+6]  = 1.
	Memr[xcoeff+7]  = 336.
	Memr[xcoeff+8]  =  0.4259691282581062 + 1089.395693397732 - 42. - 0.7 - 18.3		# -14.4 for mirror realignment
	Memr[xcoeff+9]  =  1.135406932277368  + 694.8209734223686
	Memr[xcoeff+10] =  0.4332245141494603
	Memr[xcoeff+11] =  0.6601396941753989
	Memr[xcoeff+12] = -0.364762294701932 + 2.9168802655806751
	Memr[xcoeff+13] = -1.114474752798539
	Memr[xcoeff+14] = -0.07100544020606963
	Memr[xcoeff+15] =  0.01848460943466597
	Memr[xcoeff+16] =  1.418819558933747
	Memr[xcoeff+17] =  3.289284015335439
	Memr[xcoeff+18] =  0.08722518142092386
	Memr[xcoeff+19] = -0.06158022096394576

	Memr[ycoeff  ]  = 2.
	Memr[ycoeff+1]  = 3.
	Memr[ycoeff+2]  = 4.
	Memr[ycoeff+3]  = 1.
	Memr[ycoeff+4]  = 1.
	Memr[ycoeff+5]  = 216.
	Memr[ycoeff+6]  = 1.
	Memr[ycoeff+7]  = 336.
	Memr[ycoeff+8]  =  0.464683175237381  + 1096.128797371226  - 6.0 - 70.3			# -67.2 for mirror realignment
	Memr[ycoeff+9]  =  0.1765672348323415 +  1.460280443638984
	Memr[ycoeff+10] =  0.1613619533305995
	Memr[ycoeff+11] = -2.080077948000524 + -1093.182286712313
	Memr[ycoeff+12] = -3.252201997220546
	Memr[ycoeff+13] = -1.548886149883428
	Memr[ycoeff+14] =  1.719400498669824
	Memr[ycoeff+15] = -0.02326077701339106
	Memr[ycoeff+16] =  0.3049788335944159
	Memr[ycoeff+17] = -2.430251103657301
	Memr[ycoeff+18] = -0.4204804794749471
	Memr[ycoeff+19] =  0.2496466944483961

	call gsrestore (sfx, Memr[xcoeff], ncoeff)
	call gsrestore (sfy, Memr[ycoeff], ncoeff)

	call mfree (xcoeff, TY_REAL)
	call mfree (ycoeff, TY_REAL)
end

procedure	rad_refract (ra, dec, h, lat, pres, temp, wave, ara, adec)
#
#     Calculate the apparent declination and right ascention corrected
#     for atmospheric refraction, given the true RA and DEC, etc.
#     RA and DEC are in radians.
#
#     Parameters -    (">" input, "<" output)
#
#     (<) ara        (Real*8) Apparent Right Ascention in radians
#
#     (<) adec       (Real*8) Apparent Declination in radians
#  
#     (>) ra         (Real*8) True Right Ascention in radians
#
#     (>) dec        (Real*8) True Declination in radians
#
#     (>) h          (Real*8) Hour angle in radians
#  
#     (>) lat        (Real*8) Observer's latitude in radians
#
#     (>) pres       (Real*8) Local air pressure in mm Hg.
#
#     (>) temp       (Real*8) Local air temperature in Celsius.
#   
#     (>) wave       (Real*8) Wavelength in micrometers
#
#     (>) vp         (Real*8) Vapor pressure in mm Hg (not currently used.)
#                                         [ACP: because formula NOT in Allen?]
#     (<) zd         (Real*8) Zenith distance in degrees
#
#     (<) n          (Real*8) Atmospheric refractivity
#                                         [ACP: ie `R0'; approx. index(air)-1]
#     (<) r          (Real*8) Refraction coefficient
#
#     History:  CRO/CIT.  20 June 1988.  Original.
#  Converted to spp by A. Phillips, Nov 95
#
double	ara, adec, ra, dec, lat, pres, temp, wave
double	h, da, dd, num, tcorr, tanz
double	wavers, cosq, sinq, sina, sinz, cosz
double	zd, n, r
#

begin

# Altitude and Zenith Distance
	sina = sin (lat) * sin (dec) + cos (lat) * cos (dec) * cos (h)
	cosz = sina
	sinz = sqrt (1. - sina * sina)
	tanz = sinz / cosz

	if  (sinz != 0.) {
		sinq = cos (lat) * sin (h) / sinz
		cosq =  (cos (HALFPI - lat) - cosz * cos (HALFPI - dec)) / 
						(sinz * sin (HALFPI - dec))
	} else {
#  zenith distance is 0 so q doesn't matter. r=0 anyway.
		sinq = 0.
		cosq = 0.
	}

# refractive index at standard TP
	wavers = 1. / (wave * wave)
	n = 64.328 + (29498.1 / (146. - wavers)) +  (255.4 / (41. - wavers)) 
	n = 1.e-6 * n

# temperature and pressure correction
	tcorr = 1. + 0.003661 * temp
	num   = 720.88 * tcorr
	n = n * ((pres * (1. + (1.049 - 0.0157 * temp) * 1.e-6 * pres)) / num)
	
	zd = asin (sinz)
	zd = RADTODEG (zd)

# r, refraction
	r = n * 206265. * tanz

# changes in ra and dec in radians.
	da = DEGTORAD (r / 3600.) * sinq / cos (dec)
	dd = DEGTORAD (r / 3600.) * cosq

# corrected ra and dec
	ara =  ra  + da
	adec = dec + dd
end

# FP_COORDS: Get coords in Focal Plane system

procedure	fp_coord(ra, dec, ra0, dec0, cosa, sina, xx, yy, adc)

double	ra, dec				# refracted coords in radians
double	ra0, dec0			# telescope coords in radians
double	cosa, sina			# PA transform coords
double	xx, yy				# x,y coords in focal plane coords
bool	adc				# correct scale for ADC?

double	eta, nu, denom

double	delra, sindr, cosdr		# These could probably be real
double	sind, cosd, sind0, cosd0

begin
	delra = ra - ra0
	sindr = sin (delra)
	cosdr = cos (delra)
	sind = sin (dec)
	cosd = cos (dec)
	sind0 = sin (dec0)
	cosd0 = cos (dec0)
	denom = (sind * sind0 + cosd * cosd0 * cosdr)

# eta, nu in arcsec
	eta = RADTODEG(cosd * sindr / denom) * 3600.
	nu  = RADTODEG((sind * cosd0 - cosd * sind0 * cosdr) / denom) * 3600.

# delx, dely in mm wrt telescope center
	xx = ( cosa * eta + sina * nu) * MM_ARCS
	yy = (-sina * eta + cosa * nu) * MM_ARCS

# NEW! correction for ADC plate scale change:
	if (adc) {
		xx = xx * ADC_ADJ
		yy = yy * ADC_ADJ
	}
end

procedure	lris_coord (xfp, yfp, xref, yref, x, y)

double	xfp, yfp		# x,y in focal-plane system
double	xref, yref		# x,y for reference
double	x, y			# returned x,y values (wrt mask center)

double	xx, yy, yccd

begin
	xx =  (xfp - FP_XCEN) / cos (DEGTORAD(M_ANGLE)) + FP_XCEN
	yy =  yfp / cos (DEGTORAD(M_BEND))	# Crude

## This is VERY kludgy, because Judy converts to pixels to apply an
## emperically determined correction:

	x = MASK_X0 + (FP_XCEN - xx) / SCL_FIX		# These in mm, and ...
	yccd = CCDYCEN + (FP_YCEN + yy) / MM_PIX	# ... these in pix

	x = x + (ACORR1 + yccd*(ACORR2 + yccd*(ACORR3 + yccd*ACORR4))) * MM_PIX
	y = MASK_Y0 - yy - FP_YCEN

# New system
	x = x - xref
	y = y - yref
end



#
# MAP_SORT: Sort the map info, separating the file by SLIT objects, REF objects
# and the Guide Star.  This task could be a lot cleaner if we stored the
# table as an array rather than individual vectors.
#

int	procedure map_sort (y, npt, x, yb, yt, ndx, pcode, guide_star)

real	y[npt]					# y vector (to sort by)
int	npt					# length
real	x[npt]					# x vector (to sort by)
real	yb[npt], yt[npt]			# y-limit info
int	ndx[npt]				# index to ID
int	pcode[npt]				# priority code
int	guide_star				# Is there a guide star?

int	i, i2
int	code
int	ndx2
int	nref					# number of reference objects
int	nobj					# number of program targets

begin
# We start by separating the reference objects and guide star from the others:

# if guide star present, adjust so it has space at end
	if (guide_star == YES) {
		i2 = npt - 1
	} else {
		i2 = npt
	}


# now go through list backwards, switching reference,guide objects towards end
	nref = 0
	nobj = 0
	ndx2 = i2
	do i = i2, 1, -1 {
		code = pcode[i]
		if (code == CODE_GS) {		# If GS, switch to end; test new
		    	call rndxswitch (y, ndx, npt, i, npt)
			code = pcode[npt]
		}
		if (code == CODE_RF) {
# if in last position, or next also ref, don't switch
			if (i < ndx2)
				call rndxswitch (y, ndx, npt, i, ndx2)
			nref = nref + 1
			ndx2 = ndx2 - 1
		} else {
			nobj = nobj + 1
		}
	}
	ndx2 = ndx2 + 1				# first ref object


# Now sort the slits
	call ndxsortr (y, ndx, nobj)

# ... and reference objects
	if (nref > 1)
		call ndxsortr (y[ndx2], ndx[ndx2], nref)

	call re_orderr (ndx, x, x, npt, YES)
	call re_orderr (ndx, yb, yb, npt, YES)
	call re_orderr (ndx, yt, yt, npt, YES)
	call re_orderi (ndx, pcode, pcode, npt, YES)

# and return value of slits
	return (nobj)
end



# SLIT_ADJUST: lengthen/shorten slits to use maximum space
# Recall, these are based on delta values so we can treat the bar sensibly
# Only objects which fall on mask included in list, and guide star is at end
# At this point, all values in mm

procedure	slit_adjust (y, ybot, ytop, y1, y2, pcode, n, sep_slit)

real	y[n]				# y-object
real	ybot[n], ytop[n]		# min. y length, bottom and top
real	y1[n], y2[n]			# actual ymin, ymax
int	pcode[n]			# priority/type code
int	n				# length of object list
real	sep_slit			# spacing between slits

int	i
int	imax
real	ymax, del
real	wt1, wt2, sum			# weights to distribute del

begin
	ymax = y[1]
	imax = 1
	do i = 2, n {
		if (y[i] > ymax) {
			ymax = y[i]
			imax = i
		}
	}

	if (pcode[1] != CODE_AS)
		y2[1] = MIN_MY - MASK_Y0	# put at edge; ref'd to MASK_Y0
	if (pcode[imax] != CODE_AS)
		y1[imax] = MAX_MY - MASK_Y0

	do i = 2, imax {
		del = y2[i] - y1[i-1] - sep_slit

# check for alignment stars (guide star should come at end)
		if (pcode[i] == CODE_AS)
			wt2 = 0.
		else if (del < 0.)
			wt2 = 0.5
		else
			wt2 = ytop[i]

		if (pcode[i-1] == CODE_AS)
			wt1 = 0.
		else if (del < 0.)
			wt1 = 0.5
		else
			wt1 = ybot[i-1]

		sum = wt1 + wt2
		if (sum == 0.)				# both AS, no changes
			next

		y2[i]   = y2[i]   - del * wt2 / sum	# Note sign!
		y1[i-1] = y1[i-1] + del * wt1 / sum
	}

# Now adjust at the bar so that no slit crosses the mask center:
	do i = 2, imax {
		if (y[i] > 0. && y[i-1] < 0.) {
		    if (pcode[i] != CODE_AS && pcode[i] != CODE_GS)
			y2[i] = 0.
		    if (pcode[i-1] != CODE_AS && pcode[i-1] != CODE_GS)
			y1[i-1] = 0.
		}
	}
end
# Ideally, we should keep a buffer of "reserves" for each end, so that
# we can adjust for a "min_length" condition, too.

#
# RE_ORDERR: reorder a real vector by index (input and output may be same)
#
procedure re_orderr (ndx, in, out, n, zero_index)

int	ndx[n]				# index of positions
real	in[n]				# input vector to be reordered
real	out[n]				# output vector -- reodered by index
int	n				# length of vectors
int	zero_index			# Is the vector zero-indexed?

int	i
int	ioff				# offset for zero indexing
pointer	sp, buf				# pointer for stack, work vector

begin
# allocate stack memory for work vector:
	call smark (sp)
	call salloc (buf, n, TY_REAL)

# do the reordering:
	if (zero_index == YES)
		ioff = 1
	else
		ioff = 0

	do i = 1, n
		Memr[buf+i-1] = in[ndx[i]+ioff]

# copy to output vector
	call amovr (Memr[buf], out, n)

	call sfree (sp)
end

#
# RE_ORDERI: reorder an integer vector by index (input and output may be same)
#
procedure re_orderi (ndx, in, out, n, zero_index)

int	ndx[n]				# index of positions
int	in[n]				# input vector to be reordered
int	out[n]				# output vector -- reodered by index
int	n				# length of vectors
int	zero_index			# Is the vector zero-indexed?

int	i
int	ioff				# offset for zero indexing
pointer	sp, buf				# pointer for stack, work vector

begin
# allocate stack memory for work vector:
	call smark (sp)
	call salloc (buf, n, TY_INT)

# do the reordering:
	if (zero_index == YES)
		ioff = 1
	else
		ioff = 0

	do i = 1, n
		Memi[buf+i-1] = in[ndx[i]+ioff]

# copy to output vector
	call amovi (Memi[buf], out, n)

	call sfree (sp)
end

#
# NDXSORTR: Sort a list of real values (low to high) and carry index
# (see also re_orderX)

procedure	ndxsortr (r, ndx, n)

real	r[n]				# Real vector to sort
int	ndx[n]				# index vector
int	n				# number of points

int	i, j
real	h
int	ih

begin
# Sort the list (low-to-high)
	do i = 1, n-1 {
	    do j = 1, n-i {
		if (r[j] > r[j+1]) {
			h = r[j+1]
			r[j+1] = r[j]
			r[j] = h

			ih = ndx[j+1]
			ndx[j+1] = ndx[j]
			ndx[j] = ih
		}
	    }
	}
end
#
# RNDXSWITCH: Switch 2 lines of real items and index 
#

procedure	rndxswitch (r, ndx, n, i1, i2) 

real	r[n]				# Real vector
int	ndx[n]				# Index vector
int	n				# number of points
int	i1, i2				# indices of items to be switched

real	h
int	ih

begin
	h = r[i1]
	r[i1] = r[i2]
	r[i2] = h

	ih = ndx[i1]
	ndx[i1] = ndx[i2]
	ndx[i2] = ih
end

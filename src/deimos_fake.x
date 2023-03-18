include <math.h>
include <gset.h>
include <gim.h>
include "futil.h"

define		CCD_NX	2048		# X-Size of CCD (pix)
define		CCD_NY	6000		# X-Size of CCD (pix)
define		CCD_X0	1024		# X-Center of CCD (pix)
define		CCD_Y0	3000		# Y-Center of CCD (pix)
define		GUID_X1	 20		# Lower X-limit to guider box (mm)
define		GUID_X2	196		# Upper X-limit to guider box (mm)
define		GUID_Y1	386		# Lower Y-limit to guider box (mm)
define		GUID_Y2	416		# Upper Y-limit to guider box (mm)
define		ASECPMM	1.378		# Arcsec/mm
define		ASECPPX	0.215		# Arcsec/pix
define		OBS_LAT	19.8		# Keck latitude
define		CAM_ANG 44.08		# Angle of camera wrt collimator
define		GEN_ZMX	75.0		# General max zenith angle
define		SRV_AZ1	185.0		# Azimuth of Service tower exclusion K2
define		SRV_AZ2	332.0		# Azimuth of Service tower exclusion K2
define		SRV_ZMX	53.2		# Z-angle of Service tower exclusion K2
# define	SRV_AZ1	1.6		# Azimuth of Service tower exclusion
# define	SRV_AZ2	151.0		# Azimuth of Service tower exclusion
# define	SRV_ZMX	54.		# Z-angle of Service tower exclusion

define	KEYSFILE	"ucsclris$lib/keys/simulator.key"

define		MASK_X0	88.0 		# X-Center of mask (mm)
define		MASK_Y0	360.0  		# Y-Center of mask (mm)
define		XSENSE	+1.0		# X sense of mask coord system
define		YSENSE	-1.0		# Y sense of mask coord system

define		PA_INCR	1.		# PA Increment (deg) for p/n commands
define		TR_INCR	1.		# Translational Incr. (") for hjkl
define		SZ_INLN	132		# Maximum length of input line
define		SLITLEN 8.		# Slit length (") in disp. diagram
define		PRIMRY	1		# Primary object flag
define		SECDRY	2		# Secondary object flag
define		ALIGN	3		# Secondary object flag
define		SZ_ID	9		# Length of input ID field (char)

# Mode dependent:
define		MDEF_WD	880.		# Mask Default field width (arcsec)
define		LDEF_WD	250.		# Longslit Default field width (arcsec)
define		MDEF_XM	1.		# Mask Default x-magn factor
define		LDEF_XM	4.		# Longslit Default x-magn factor

# Mask specific:
#define		BAR_LOW	161.3		# Lower Bar Y (mm)
#define		BAR_UPP	169.4		# Upper Bar Y (mm)
define		BAR_LOW	164.6		# Lower Bar Y (mm)
define		BAR_UPP	171.0		# Upper Bar Y (mm)
define		MIR_X0	46.3		# Mirror offset along Bar X (mm)
define		MIR_Y0	-21.0		# Mirror (#2) offset from Bar Y (mm)
define		MIR_XSZ	22.9		# Mirror size in X (mm)
define		MIR_YSZ	22.9		# Mirror size in Y (mm)

# Longslit specific:
define		SLIT_Y1	80		# Lower Y-limit to slit (mm)
define		SLIT_Y2	240		# Upper Y-limit to slit (mm)

define		NLIMITS	34		# number of limit parameters following
define		CX1	$1[1]		# CCD Lower X-limit (" rel. center)
define		CY1	$1[2]		# CCD Lower Y-limit (" rel. center)
define		CX2	$1[3]		# CCD Upper X-limit (" rel. center)
define		CY2	$1[4]		# CCD Upper Y-limit (" rel. center)
define		MX1	$1[5]		# Mask Lower X-limit (" rel. center)
define		MY1	$1[6]		# Mask Lower Y-limit (" rel. center)
define		MX2	$1[7]		# Mask Upper X-limit (" rel. center)
define		MY2	$1[8]		# Mask Upper Y-limit (" rel. center)
define		BY1	$1[9]		# Bar Lower Y (" rel. center)
define		BY2	$1[10]		# Bar Upper Y (" rel. center)
define		GX1	$1[11]		# Guider Lower X-limit (" rel. center)
define		GY1	$1[12]		# Guider Lower Y-limit (" rel. center)
define		GX2	$1[13]		# Guider Upper X-limit (" rel. center)
define		GY2	$1[14]		# Guider Upper Y-limit (" rel. center)
define		X1	$1[15]		# Eff. Lower X-limit (" rel. center)
define		Y1	$1[16]		# Eff. Lower Y-limit (" rel. center)
define		X2	$1[17]		# Eff. Upper X-limit (" rel. center)
define		Y2	$1[18]		# Eff. Upper Y-limit (" rel. center)
define		DISP1	$1[19]		# Starting dispersion (")
define		DISP2	$1[20]		# Ending dispersion (")
define		PA1	$1[21]		# Starting Disp. PA (deg)
define		PA2	$1[22]		# Ending Disp. PA (deg)
define		ELEV_1	$1[23]		# Starting Dewar Elevation (deg)
define		ELEV_2	$1[24]		# Ending Dewar Elevation (deg)
define		SLITWID	$1[25]		# Width of slit (")
define		MINSLIT	$1[26]		# Slit length in plot (")
define		ERR_LOW	$1[27]		# If non-zero, service tower conflict
define		ERR_LN2	$1[28]		# If non-zero, LN2 spillage
define		PL_XMAG	$1[29]		# Xmag in plot (default 1)
define		PL_FWID	$1[30]		# Field width (")
define		MIR_X1	$1[31]		# Mirror Lower X-limit (" rel. center)
define		MIR_X2	$1[32]		# Mirror Upper X-limit (" rel. center)
define		MIR_Y1	$1[33]		# Mirror Lower Y-limit (" rel. center)
define		MIR_Y2	$1[34]		# Mirror Upper Y-limit (" rel. center)

# T_SIMULATOR: LRIS simulator; formerly...
#    T_MSSELECT: Multi-slit select; formerly ...
#    T_AUTOSLUT: Because it panders to AUTOSLIT (and it's user friendly, too)
#
# Note on Coord system:
# The x axis is parallel to the slit-mask-punch x-axis (ie neg fp-coord system)
# The y-axis is anti-parallel to the slit-punch y-axis (ie pos fp-coord system)
# This produces an image as seen on the CCD
# All coords relevant to the mask (mask, bar, pickoff mirrors) are given in
# slit-punch coords (mm).
# All displayed coords are in arc-sec, with a Mercator-style dec value
# (normalized to the central dec value).

procedure t_deimos_fake()

# Updates
# -- New format (slit top/bottom)
# -- CENTER line in input
# -- Alignment/Guide star select

# These all need to be incorporated:
# -- Spectral range			# (done) get accurate CCD pos.
# -- Bar				# Done
# -- Checks on PA (for LN2 spillage)	# Done
# -- Check on dispersion angle		# Done
# -- Check on tower exclusion		# Done
# -- Punch machine x-limit		# done
# -- Guide star possibilities		# (done) get box location

# Current approximations:
# -- edges are straight (ie won't work at pole)
# -- All coord. epochs are current

# Needs to have the following options:
# -- Change field center		# Y
# -- Change PA				# Y
# -- Add/delete objects			# (Y)
# -- provide info on particular objects # (Y)
# -- set size of viewing region		# Y

bool	long_slit_mode
char	objfile[SZ_FNAME]			# ID, prior, mag, RA, Dec
char	other[SZ_FNAME]				# ID, prior, mag, RA, Dec
char	output[SZ_FNAME]			# Output file name
char	remain[SZ_FNAME]			# Output list of remaining obj.
double	dra0, ddec0				# Coords. of field center
double	equinox					# Equinox of Coords
real	limit[NLIMITS]
pointer	fda, fdb, fdc, fdd

char	tchar
double	dra, ddec
double	cosd, cosd0
char	id[SZ_ID]
int	priority
real	magn, center
int	i
int	npt, ndx
real	x1, x2, y1, y2
real	wave1, wave2, disp, woff
real	ha1, ha2
real	dec, lat, pa, atmdisp, slit, z1, z2, az1, az2
real	xoff, yoff, theta
pointer	bufx, bufy, bufxt, bufyt, buff, bufin, bufp
real	slitbot, slittop
pointer	bufpa, bufs1, bufs2		# PA, slit_below, slit_above

real	chk_elev()

bool	clgetb(), strne(), streq()
double	clgetd()
real	clgetr()
int	fscan(), nscan()
pointer	open()

begin
	call clgstr ("objfile", objfile, SZ_FNAME)
	call clgstr ("other_obj", other, SZ_FNAME)
	call clgstr ("output", output, SZ_FNAME)
	call clgstr ("remain", remain, SZ_FNAME)
	ha1 = clgetr ("ha")
	ha2 = clgetr ("exposure") + ha1			# Ending HA
	slit = clgetr ("slit")
	wave1 = clgetr ("blue")
	wave2 = clgetr ("red")
	disp = clgetr ("dispersion")
	long_slit_mode = clgetb ("long_slit_mode")

	if (!long_slit_mode) {
		x1 = clgetr ("x1")
		x2 = clgetr ("x2")
		y1 = clgetr ("y1")
		y2 = clgetr ("y2")
		woff = 0.5 * abs (wave2-wave1) / disp
	} else {
		woff = 0.
	}

# Open the primary file; check for CENTER specifications
        fda = open (objfile, READ_ONLY, TEXT_FILE)
	while (fscan(fda) != EOF) {
		call gargwrd (tchar, 1)
		if (tchar == '#' || nscan() == 0) {
			next
		}
		call reset_scan()

		call gargwrd (id, SZ_ID)
		if (streq (id, "CENTER")) {
#	if present, get info
			call gargi (priority)
			call gargd (equinox)
			call gargd (dra0)
			call gargd (ddec0)
			call gargr (theta)
			if (nscan() < 6 || theta == INDEF)
				theta = 0.
			if (nscan() >= 5) {
			    call printf ("Field Center, Eqx. from Input File\n")
			    break
			} else {
			    call eprintf ("Incorrect CENTER specification\n")
			}
		}
#	otherwise get from parameters (forced)
		dra0 = clgetd ("ra0")
		ddec0 = clgetd ("dec0")
		equinox = clgetd ("equinox")
		theta = clgetr ("PA")
		break
	}
	call seek (fda, BOF)

# Set the limits (mostly from definitions above), in image coord space
# CCD:
	CX1(limit) = (1 - CCD_X0 + woff) * ASECPPX
	CX2(limit) = (CCD_NX - CCD_X0 - woff) * ASECPPX
	CY1(limit) = (1 - CCD_Y0) * ASECPPX 
	CY2(limit) = (CCD_NY - CCD_Y0) * ASECPPX 
# Guider:
	GX1(limit) = (GUID_X1 - MASK_X0) * ASECPMM
	GX2(limit) = (GUID_X2 - MASK_X0) * ASECPMM
	GY1(limit) = (GUID_Y1 - MASK_Y0) * ASECPMM 
	GY2(limit) = (GUID_Y2 - MASK_Y0) * ASECPMM 

	if (long_slit_mode) {
# Longslit (treat as mask):
		MX1(limit) = -0.5 * slit
		MX2(limit) =  0.5 * slit
		if (YSENSE > 0.) {
			MY1(limit) = (SLIT_Y1 - MASK_Y0) * ASECPMM 
			MY2(limit) = (SLIT_Y2 - MASK_Y0) * ASECPMM 
		} else {
			MY2(limit) = -(SLIT_Y1 - MASK_Y0) * ASECPMM 
			MY1(limit) = -(SLIT_Y2 - MASK_Y0) * ASECPMM 
		}
# Bar and pick-off mirror (not there!)
		BY1(limit) =  GY1 (limit)
		BY2(limit) =  GY1 (limit)
		MIR_X1(limit) = GX1 (limit)
		MIR_X2(limit) = GX1 (limit)
		MIR_Y1(limit) = GY1 (limit)
		MIR_Y2(limit) = GY1 (limit)
	} else {
# Slit mask:
		if (XSENSE > 0.) {
			MX1(limit) = (x1 - MASK_X0) * ASECPMM
			MX2(limit) = (x2 - MASK_X0) * ASECPMM
		} else {
			MX2(limit) = -(x1 - MASK_X0) * ASECPMM
			MX1(limit) = -(x2 - MASK_X0) * ASECPMM
		}
		if (YSENSE > 0.) {
			MY1(limit) = (y1 - MASK_Y0) * ASECPMM 
			MY2(limit) = (y2 - MASK_Y0) * ASECPMM 
			BY1(limit) = (BAR_LOW - MASK_Y0) * ASECPMM 
			BY2(limit) = (BAR_UPP - MASK_Y0) * ASECPMM 
		} else {
			MY2(limit) = -(y1 - MASK_Y0) * ASECPMM 
			MY1(limit) = -(y2 - MASK_Y0) * ASECPMM 
			BY2(limit) = -(BAR_LOW - MASK_Y0) * ASECPMM 
			BY1(limit) = -(BAR_UPP - MASK_Y0) * ASECPMM 
		}
# Slit-viewing Pickoff mirror location:
		center = XSENSE * MIR_X0
		MIR_X1(limit) = (center - 0.5*MIR_XSZ) * ASECPMM
		MIR_X2(limit) = MIR_X1(limit) + MIR_XSZ * ASECPMM
		center = (BY1(limit) + BY2(limit)) * 0.5 + MIR_Y0
		MIR_Y1(limit) = (center - 0.5*MIR_YSZ) * ASECPMM
		MIR_Y2(limit) = MIR_Y1(limit) + MIR_YSZ * ASECPMM
	}

# Effective slit limits:
	X1(limit) = max (CX1(limit), MX1(limit))
	Y1(limit) = max (CY1(limit), MY1(limit))
	X2(limit) = min (CX2(limit), MX2(limit))
	Y2(limit) = min (CY2(limit), MY2(limit))

# Dispersion info:
	dec = real (ddec0)
	lat = OBS_LAT
	call atm_geom (ha1, dec, wave1, wave2, lat, z1, az1, pa, atmdisp)
	PA1(limit) = pa
	DISP1(limit) = atmdisp
	call atm_geom (ha2, dec, wave1, wave2, lat, z2, az2, pa, atmdisp)
	PA2(limit) = pa
	DISP2(limit) = atmdisp
	if (dec > lat) {
	    if (az1 > 180.)
		az1 = az1 - 360.
	    if (az2 > 180.)
		az2 = az2 - 360.
	} else {
	    if (az1 < 0.) 
		az1 = az1 + 360.
	    if (az2 < 0.) 
		az2 = az2 + 360.
	}
	SLITWID(limit) = slit

# Check limits:
	ERR_LOW(limit) = chk_elev (az1, z1, az2, z2)
	ERR_LN2(limit) = 0.
	if (ERR_LOW(limit) != 0.)
		call printf ("Potential Service Tower/Shutter conflict!\n")

# Initial plot limits:
	if (long_slit_mode) {
		PL_FWID(limit) = LDEF_WD
		PL_XMAG(limit) = LDEF_XM
	} else {
		PL_FWID(limit) = MDEF_WD
		PL_XMAG(limit) = MDEF_XM
	}
	MINSLIT(limit) = clgetr ("min_slit")

# Count the primary list (already open) and the secondary list if present
	npt = 0
	while (fscan(fda) != EOF)
		npt = npt + 1
	call seek (fda, BOF)

	if (strne (other, "")) {
        	fdb = open (other, READ_ONLY, TEXT_FILE)
		while (fscan(fdb) != EOF)
			npt = npt + 1
		call seek (fdb, BOF)
	}

# Allocate memory
	call malloc (bufx, npt, TY_REAL)
	call malloc (bufy, npt, TY_REAL)
	call malloc (bufxt, npt, TY_REAL)
	call malloc (bufyt, npt, TY_REAL)
	call malloc (bufp, npt, TY_INT)			# priority
	call malloc (bufpa, npt, TY_REAL)		# PA
	call malloc (bufs1, npt, TY_REAL)		# slit len below
	call malloc (bufs2, npt, TY_REAL)		# slit len above
	call malloc (buff, npt, TY_INT)			# flag
	call malloc (bufin, npt*SZ_INLN, TY_CHAR)
	call amovkc (EOS, Memc[bufin], npt*SZ_INLN)

# Get the primary list (avoid center if present)
	ndx = 0
	cosd0 = cos (DEGTORAD (ddec0))
	while (fscan(fda) != EOF) {
		call gargwrd (tchar, 1)
		if (tchar == '#' || nscan() == 0) {
			next
		}
		call reset_scan()
		call gargwrd (id, SZ_ID)
		if (ndx == 0 && streq (id, "CENTER"))
			next
		call gargi (priority)
		call gargr (magn)
		call gargd (dra)
		call gargd (ddec)

		if (nscan() < 5) {
		    call eprintf ("Poorly-formatted data line %d -- skipped\n")
			call pargi (ndx+1)
		    next
		}
		cosd = cos (DEGTORAD (ddec))
		Memr[bufx+ndx] = (ddec - ddec0) * 3600. * cosd0 / cosd
		Memr[bufy+ndx] = (dra - dra0) * 15. * 3600. * cosd
		Memi[bufp+ndx] = priority

		call gargr (pa)
		call gargr (slitbot)
		call gargr (slittop)

		if (nscan() < 6) {
			pa = INDEF
		}
		if (nscan() < 7 || slitbot <= 0. || slitbot == INDEF)
			slitbot = MINSLIT(limit) * 0.5
		if (nscan() < 8 || slittop <= 0. || slittop == INDEF)
			slittop = MINSLIT(limit) * 0.5
		Memr[bufpa+ndx] = pa
		Memr[bufs1+ndx] = slitbot
		Memr[bufs2+ndx] = slittop
		
		call reset_scan()
		call gargstr (Memc[bufin+ndx*SZ_INLN], SZ_INLN)
		if (priority >= 0)
			Memi[buff+ndx] = -1 * PRIMRY
		else
			Memi[buff+ndx] = -1 * ALIGN
		ndx = ndx + 1
	}
	call close (fda)

# Get the secondary list:
	if (strne (other, "")) {
	    while (fscan(fdb) != EOF) {
		call gargwrd (tchar, 1)
		if (tchar == '#' || nscan() == 0) {
			next
		}
		call reset_scan()
		call gargwrd (id, SZ_ID)
		call gargi (priority)
		call gargr (magn)
		call gargd (dra)
		call gargd (ddec)

		if (nscan() < 5) {
		    call eprintf ("Poorly-formatted data line %d -- skipped\n")
			call pargi (ndx+1)
		    next
		}
		cosd = cos (DEGTORAD (ddec))
		Memr[bufx+ndx] = (ddec - ddec0) * 3600. * cosd0 / cosd
		Memr[bufy+ndx] = (dra - dra0) * 15. * 3600. * cosd
		Memi[bufp+ndx] = priority

		call gargr (pa)
		call gargr (slitbot)
		call gargr (slittop)

		if (nscan() < 6) {
			pa = INDEF
		}
		if (nscan() < 7 || slitbot <= 0. || slitbot == INDEF)
			slitbot = MINSLIT(limit) * 0.5
		if (nscan() < 8 || slittop <= 0. || slittop == INDEF)
			slittop = MINSLIT(limit) * 0.5
		Memr[bufpa+ndx] = pa
		Memr[bufs1+ndx] = slitbot
		Memr[bufs2+ndx] = slittop
		
		call reset_scan()
		call gargstr (Memc[bufin+ndx*SZ_INLN], SZ_INLN)
		if (priority >= 0)
			Memi[buff+ndx] = -1 * SECDRY
		else
			Memi[buff+ndx] = -1 * ALIGN
		ndx = ndx + 1
	    }
	    call close (fdb)
	}

	npt = ndx

# Open the output file
	if (strne (output, ""))
        	fdc = open (output, NEW_FILE, TEXT_FILE)
	if (strne (remain, ""))
        	fdd = open (remain, NEW_FILE, TEXT_FILE)

	call focpl_graph (Memr[bufx], Memr[bufy], Memr[bufxt], Memr[bufyt],
		Memr[bufpa], Memr[bufs1], Memr[bufs2], Memi[bufp], Memi[buff],
		npt, limit, dec, ha1, ha2, theta, xoff, yoff, bufin)

	ddec = ddec0 + xoff/3600.
	cosd = cos (DEGTORAD (ddec))
	dra = dra0 + yoff/3600./15. / cosd

	call printf ("\n\nField center: %11.1h %11.0h,  PA: %5.1f\n\n")
		call pargd (dra)
		call pargd (ddec)
		call pargr (theta)
	call printf ("Geometry: z1 = %4.1f (az%4.0f), z2 = %4.1f (az%4.0f)\n")
		call pargr (z1)
		call pargr (az1)
		call pargr (z2)
		call pargr (az2)
	call printf ("Dispersion: %5.2f'' @%6.1f deg --%5.2f'' @%6.1f deg\n\n")
		call pargr (DISP1(limit))
		call pargr (PA1(limit))
		call pargr (DISP2(limit))
		call pargr (PA2(limit))
	call printf ("Dewar elevation: %5.1f (start)  %5.1f (end)\n")
		call pargr (ELEV_1(limit))
		call pargr (ELEV_2(limit))

	if (ERR_LOW(limit) == 1.)
		call printf ("PROBLEM: Service Tower Conflict\n")
	if (ERR_LOW(limit) == 2.)
		call printf ("PROBLEM: Occulted by Shutter\n")
	if (ERR_LN2(limit) != 0.)
		call printf ("PROBLEM: PA May Cause LN2 Spillage\n")

# If output file specified, copy same info.
	if (strne (output, "")) {
	    call fprintf (fdc,
		"#\n#\n# Field center: %11.1h %11.0h,  PA: %7.3f\n#\n")
		call pargd (dra)
		call pargd (ddec)
		call pargr (theta)
	    call fprintf (fdc,
		"# Geometry: z1 = %4.1f (az%4.0f), z2 = %4.1f (az%4.0f)\n")
		call pargr (z1)
		call pargr (az1)
		call pargr (z2)
		call pargr (az2)
	    call fprintf (fdc,
		"# Dispersion: %5.2f'' @%6.1f deg --%5.2f'' @%6.1f deg\n#\n")
		call pargr (DISP1(limit))
		call pargr (mod (PA1(limit), 180.))
		call pargr (DISP2(limit))
		call pargr (mod (PA2(limit), 180.))
	    call fprintf (fdc,
			"# Dewar elevation: %5.1f (start)  %5.1f (end)\n")
		call pargr (ELEV_1(limit))
		call pargr (ELEV_2(limit))

	    if (ERR_LOW(limit) == 1.)
		call fprintf (fdc, "# PROBLEM: Service Tower Conflict\n")
	    if (ERR_LOW(limit) == 2.)
		call fprintf (fdc, "# PROBLEM: Occulted by Shutter\n")
	    if (ERR_LN2(limit) != 0.)
		call fprintf (fdc, "# PROBLEM: PA May Cause LN2 Spillage\n")
	}
	
# Output list to file or STDOUT:
	if (strne (output, "")) {
	    call fprintf (fdc, "\nCENTER   9999 %7.2f %11.2h %11.1h %5.1f\n")
		call pargd (equinox)
		call pargd (dra)
		call pargd (ddec)
		call pargr (theta)
	    do i = 0, npt-1 {
		if (Memi[buff+i] > 0) {
		    call fprintf (fdc, "%s\n")
			call pargstr (Memc[bufin+i*SZ_INLN])
		}
	    }
	    call close (fdc)
	} else {
	    call printf ("\nCENTER   9999 %7.2f %11.2h %11.1h %5.1f\n")
		call pargd (equinox)
		call pargd (dra)
		call pargd (ddec)
		call pargr (theta)
	    do i = 0, npt-1 {
		if (Memi[buff+i] > 0) {
		    call printf ("%s\n")
			call pargstr (Memc[bufin+i*SZ_INLN])
		}
	    }
	}

# If remainder file specified, write out unused objects
	if (strne (remain, "")) {
	    do i = 0, npt-1 {
		if (Memi[buff+i] < 0) {
		    call fprintf (fdd, "%s\n")
			call pargstr (Memc[bufin+i*SZ_INLN])
		}
	    }
	    call close (fdd)
	}

	call mfree (bufin, TY_CHAR)
	call mfree (bufp, TY_INT)
	call mfree (buff, TY_INT)
	call mfree (bufs2, TY_REAL)
	call mfree (bufs1, TY_REAL)
	call mfree (bufpa, TY_REAL)
	call mfree (bufyt, TY_REAL)
	call mfree (bufxt, TY_REAL)
	call mfree (bufy, TY_REAL)
	call mfree (bufx, TY_REAL)
end


#####
## TBD:
# 1. all logic; probably need to parse window, IMTYPE
# 2. Make sure keywords are present
# 3. Add binning
# 4. Add windowing? (Currently assumes NX,NY are full detector in amp_parse)
# 5. Add scaling
#####
##### Blue
## WINDOW  = '1,0,0,2048,4096'
## PREPIX  =                   51
## POSTPIX =                   80
## IMTYPE  = 'n MOSAIC'
## BINNING = '1,1     '
## AMPPSIZE= '[1:1024,1:4096]'
## DETLSIZE= '[1:4096,1:4096]'
######
#####cf RED:
## AMPLIST = '2,1,0,0 '
## NUMAMPS =                    2
## REDXFLIP=                    T
## IMTYPE  = 'TWOAMPTOP'
## PREPIX  =                    0
## PRELINE =                    0
## POSTPIX =                   10
## YFLIP   =                    1
## BINNING = '1, 1    '
## CCDPSIZE= '[1:2048,1:2048]'
## WINDOW  = '0,0,0,2048,2048'
############

include	<imhdr.h>
include	<error.h>

define	RED	1	# red side
define	BLUE	2	# blue side
define	MAX_AMP	4	# Maximum number of amplifiers

define	SZ_KEYW	8	# character size of keyword
define	SZ_KVAL	80	# character size of keyword value
 
# T_L4PROCESS: Do a 1-pass bias/trim/flat-field, allowing for up to 4 amplifiers
#
#  Needs a check for bwidth, cskip inside bounds! (from l2proc)

define		SCALE	0.5		# Scaling for use with SHORT format

procedure t_l4process()

char	image[SZ_FNAME]
char	flatim[SZ_FNAME]
char	imageout[SZ_FNAME], imorig[SZ_FNAME]
char	history[SZ_LINE]
bool	flatten					# flatten the image?
bool	intscale				# half and store as short?
#int	x1, y1, x2, y2				# region of interest
int	l1, l2					# lines of bias strip
#int	low1, low2				# col. of lower bias strip
#int	upp1, upp2				# col. of upper bias strip
#int	nxlow					# last x of lower half of CCD
pointer	im1, im2, im3
long	v1[IM_MAXDIM], v2[IM_MAXDIM], v3[IM_MAXDIM]
pointer	line1, line2, line3

int	nx, ny					# size of output image
# real	bias1, bias2				# bias levels 1 and 2
# int	xoff, i
int	i, j
int	bwid, bmid
real	bsum
real	scaling
pointer	bufb

# new
char	side[SZ_LINE]		# side (red/blue)
int	namp			# number of amps
int	acol1[MAX_AMP]		# amplifier image start col
int	bcol1[MAX_AMP]		# amplifier bias  start col
int	awid			# width of active, bias regions
real	bias[MAX_AMP]		# bias levels
int	ndx

bool	strne()
bool	clgetb(), streq()
int	clgeti(), impnlr(), imgnlr(), asoki()
pointer	immap(), imgs2i()

begin
	call clgstr ("image", image, SZ_FNAME)
	call clgstr ("out_image", imageout, SZ_FNAME)
	call clgstr ("flatim", flatim, SZ_FNAME)
	flatten = strne (flatim, "")
	intscale = clgetb ("intscale")

#	x1 = clgeti ("x1")
#	y1 = clgeti ("y1")
#	x2 = clgeti ("x2")
#	y2 = clgeti ("y2")
#	nxlow = clgeti ("nxlow")
	l1 = clgeti ("l1")
	l2 = clgeti ("l2")
#	low1 = clgeti ("low1")
#	low2 = clgeti ("low2")
#	upp1 = clgeti ("upp1")
#	upp2 = clgeti ("upp2")

	im1 = immap (image, READ_ONLY, 0)
	
	call clgstr ("side", side, SZ_LINE)
call eprintf ("Side=%s\n"); call pargstr (side)
	if (streq (side, "blue")) {
		call amp_parse (im1, BLUE, nx, ny, namp, acol1, bcol1, awid, bwid)
	} else if (streq (side, "red")){
		call amp_parse (im1, RED,  nx, ny, namp, acol1, bcol1, awid, bwid)
	} else {
		call fatal (0, "Unknown side X")
	}

# Get columns for data, bias:

	if (mod (bwid, 2) != 1) {
		call eprintf ("WARNING: bias width decreased by 1\n")
		bwid = bwid - 1
	}
	bmid = (bwid + 1) / 2

# calculate the bias levels:
	call amovkr (INDEF, bias, MAX_AMP)
	do i = 1, namp {
		bufb = imgs2i (im1, bcol1[i], bcol1[i]+bwid-1, l1, l2)
		bsum = 0.
		do j = 0, l2-l1 {
			bsum = bsum + asoki (Memi[bufb+j*bwid], bwid, bmid)
		}
		bias[i] = bsum / (l2 - l1 + 1)
call eprintf ("bias[%d] = %6.2f\n"); call pargi (i); call pargr (bias[i])
	}

	if (flatten) {
		im2 = immap (flatim, READ_ONLY, 0)

		if (nx != IM_LEN(im2,1) || ny != IM_LEN(im2,2))
			call fatal ("Incongruous image sizes!\n")
	}

# open output image and write
	call xt_mkimtemp (image, imageout, imorig, SZ_FNAME) # ref pkg$xtools
	im3 = immap (imageout, NEW_COPY, im1)
	IM_LEN(im3,1) = nx
	IM_LEN(im3,2) = ny
	if (intscale) {
		IM_PIXTYPE(im3) = TY_SHORT
		scaling = SCALE
	} else {
		IM_PIXTYPE(im3) = TY_REAL
		scaling = 1.
	}

# add info to image header
	call sprintf (history, SZ_LINE, ": lproc (%d) bias %6.1f,%6.1f,%6.1f,%6.1f; flat=%s")
		call pargi (namp)
		do i = 1, MAX_AMP
			call pargr (bias[i])
		if (flatten)
			call pargstr (flatim)
		else
			call pargstr ("NO FLAT")
	if (intscale)
		call strcat ("; SCALED\n", history, SZ_LINE)
	else
		call strcat ("\n", history, SZ_LINE)
	call bksp_strcat (history, IM_HISTORY(im3), SZ_IMHIST)

# do the copy
	call amovkl (long(1), v1, IM_MAXDIM)
	call amovkl (long(1), v2, IM_MAXDIM)
	call amovkl (long(1), v3, IM_MAXDIM)
	v1[2] = 1	# was y1 -- WORK for sub-full frames!

	while (impnlr(im3, line3, v3) != EOF ) {
		if (imgnlr (im1, line1, v1) == EOF)
			call eprintf (" READ ERROR (input)\n")
		do i = 1, namp {
			ndx = (i-1)*awid
			call aaddkr (Memr[line1+acol1[i]-1], -bias[i], Memr[line3+ndx], awid)
		}
		if (flatten) {
		    if (imgnlr (im2, line2, v2) == EOF)
			call eprintf (" READ ERROR (flat)\n")
		    call adivr (Memr[line3], Memr[line2], Memr[line3], nx)
		} else if (scaling != 1.) {
## scaling in here
		    call amulkr (Memr[line3], scaling, Memr[line3], nx)
		}
		
	}

	if (flatten)
		call imunmap (im2)
	call imunmap (im1)
	call imunmap (im3)
	call xt_delimtemp (imageout, imorig)
end

procedure	amp_parse (im, side, nx, ny, namp, acol1, bcol1, awid, bwid)

pointer	im			# image pointer
int	side			# side (RED/BLUE)
int	nx, ny			# size of processed image
int	namp			# number of amps
int	acol1[ARB]		# array of starting columns for processing
int	bcol1[ARB]		# array of starting columns for bias regions
int	awid			# width (col) of active regions
int	bwid			# width (col) of bias   regions

int	i, ndx
int	nxdet			# full detector size
int	i1, j1, j2, stat
int	prepix, postpix

int	imgeti(), get_datasec()

begin
	if (side == BLUE) {
		stat = get_datasec (im, "DETLSIZE", i1, nxdet, j1, j2)
		stat = get_datasec (im, "AMPPSIZE", i1, awid, j1, j2)
		namp = nxdet / awid
		nx = nxdet
		ny = j2
call eprintf ("  %d amps found\n"); call pargi (namp)
	} else if (side == RED) {
		stat = get_datasec (im, "CCDPSIZE", i1, nxdet, j1, j2)
call eprintf ("RED SIDE UNFINISHED -- assumes 2 amps!!\n")
		namp = 2
		nx = nxdet
		ny = j2
		awid = nx / 2
	} else {
		call fatal (0, "Unknown LRIS side!")
	}

	prepix = imgeti (im, "PREPIX")
	postpix = imgeti (im, "POSTPIX")
	bwid = postpix

	ndx = namp * prepix + 1
	do i = 1, namp {
		acol1[i] = ndx + (i-1) * awid
		bcol1[i] = ndx + nxdet + (i-1) * postpix
	}

	if (bcol1[namp]-1+bwid != IM_LEN(im,1)) {
		call fatal (0, "KEYWORDS inconsistent with image size!")
	}
end

#
# GET_DATASEC: decode datasec.  Currently assumes NAXIS=2 (STOLEN FROM DMOS)
# In this modified form, adds keyword name to input values
#

define	SZ_FITS_STR	72

int	procedure get_datasec (im, kwnam, i1, i2, j1, j2)

pointer	im		# image descriptor
char	kwnam[ARB]	# keyword name
int	i1, i2		# datasec limits in x
int	j1, j2		# datasec limits in y

char	kwval[SZ_FITS_STR]
char	wkstr[SZ_FITS_STR]
char	tchar
int	i, n, ia, ib, io
int	naxis

int	stridx(), nscan()
begin
	naxis = IM_NDIM(im)
	if (naxis != 2)
		call fatal (0, "GET_DATASEC: naxis != 2!")

	call imgstr (im, kwnam, kwval, SZ_FITS_STR)

	tchar = '['
	ia = stridx (tchar, kwval)
	tchar = ']'
	ib = stridx (tchar, kwval)

	if (ia < 1 || ib < ia)
		call fatal (0, "GET_DATASEC: no bracket pair")

	tchar = ' '
	call amovkc (tchar, wkstr, SZ_FITS_STR)
	i = ia + 1
	io = 1
	tchar = ':'
	n = stridx (tchar, kwval[i]) - 1
	if (n < 1)
		call fatal (0, "GET_DATASEC: poor datasec format")
	call amovc (kwval[i], wkstr[io], n)

	i = i + n + 1
	io = io + n + 1
	tchar = ','
	n = stridx (tchar, kwval[i]) - 1
	if (n < 1)
		call fatal (0, "GET_DATASEC: poor datasec format")
	call amovc (kwval[i], wkstr[io], n)

	i = i + n + 1
	io = io + n + 1
	tchar = ':'
	n = stridx (tchar, kwval[i]) - 1
	if (n < 1)
		call fatal (0, "GET_DATASEC: poor datasec format")
	call amovc (kwval[i], wkstr[io], n)

	i = i + n + 1
	io = io + n + 1
	tchar = ']'
	n = stridx (tchar, kwval[i]) - 1
	if (n < 1)
		call fatal (0, "GET_DATASEC: poor datasec format")
	call amovc (kwval[i], wkstr[io], n)

	io = io + n
	wkstr[io] = EOS

	call sscan (wkstr)
		call gargi (i1)
		call gargi (i2)
		call gargi (j1)
		call gargi (j2)
	if (nscan() < 4)
		call fatal (0, "GET_DATASEC: decode failed; bad format?")
	else
		return (OK)
end

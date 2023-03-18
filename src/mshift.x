include	<math.h>
include	"lris.h"
#
# T_MSHIFT: from (x,y) (xref,yref) and PA calcuate offsets in arcsec

procedure t_mshift ()

real	x, y					# measured x and y (pixels)
real	xref, yref				# desired x and y (pixels)
real	pa					# pa in degrees
real	narcsec, earcsec			# returned offsets, N and E
bool	invert					# invert orient. (x-flipped)?
bool	dcs					# print DCS commands

char	cmdline[SZ_LINE]			# command string for DCS
int	stat					# command stat

bool	clgetb()
int	oscmd()
real	clgetr()

begin
	x = clgetr ("xobs")
	y = clgetr ("yobs")
	xref = clgetr ("xref")
	yref = clgetr ("yref")
	pa = clgetr ("pa")
	invert = clgetb ("invert")

	call mshift (x, y, xref, yref, pa, ASECPPX, narcsec, earcsec, invert)

	call printf ("====================================================\n\n")
	call printf (
		" To shift (%5.1f,%5.1f) to (%5.1f,%5.1f) at mask PA=%5.1f:\n")
		call pargr (x)
		call pargr (y)
		call pargr (xref)
		call pargr (yref)
		call pargr (pa)
	if (invert)
		call printf (" (one-amp, x-inverted mode)\n")
	call printf ("\n   MOVE TELECSOPE   %5.2f''E  and %5.2f''N \n\n")
		call pargr (earcsec)
		call pargr (narcsec)
	call printf ("====================================================\n\n")

	dcs = clgetb ("dcs")
	if (dcs) {
		call sprintf (cmdline, SZ_LINE,
		  "%s modify -s dcs2 RAOFF=%.2f DECOFF=%.2f REL2CURR=1")
			call pargstr ("rsh manuka")
			call pargr (earcsec)
			call pargr (narcsec)
		call eprintf ("\n %s \n")
			call pargstr (cmdline)
		stat = oscmd (cmdline)
		if (stat != OK) {
			call eprintf ("command failed!  (%d)\n")
				call pargi (stat)
		}

		call sprintf (cmdline, SZ_LINE, "%s waitfor -s dcs2 AXESTAT=64")
			call pargstr ("rsh manuka")	# (64=tracking)
		call eprintf ("\n %s \n")
			call pargstr (cmdline)
		stat = oscmd (cmdline)
		if (stat != OK) {
			call eprintf ("command failed!  (%d)\n")
				call pargi (stat)
		}

		call eprintf ("... done! \n")
	}
end

#
# MSHIFT: work out pixel shifts for a given MASK PA
#

procedure mshift (x, y, xref, yref, pa, pixscale, narcsec, earcsec, invert)

real	x, y					# measured x and y (pixels)
real	xref, yref				# desired x and y (pixels)
real	pa					# pa in degrees
real	pixscale				# arcsec/pixel
real	narcsec, earcsec			# returned offsets, N and E
bool	invert					# invert orient. (x-flipped)?

real	cosa, sina

begin
	cosa = cos (DEGTORAD(90.-pa))
	sina = sin (DEGTORAD(90.-pa))
	if (!invert) {
		narcsec = pixscale * ( cosa * (x-xref) + sina * (y-yref))
		earcsec = pixscale * (-sina * (x-xref) + cosa * (y-yref))
	} else {
		narcsec = pixscale * (-cosa * (x-xref) + sina * (y-yref))
		earcsec = pixscale * ( sina * (x-xref) + cosa * (y-yref))
	}
end

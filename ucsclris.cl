#{ Package UCSCLRIS is Drew Phillips' lris software


print ("\n UCSCLRIS (ver.0a for IRAF 2.12) -- Unsupported software -- User assumes risk\n")

cl < "ucsclris$lib/zzsetenv.def"
package	ucsclris, bin = ucsclrisbin$

task	maskalign,
	mboxfind,
	xbox,
	mshift,
	salign,
	flex_fit,
	l4process,
	l2process = "ucsclris$src/x_ucsclris.e"

task	$prep		= "ucsclris$prep.cl"		# Mask prep package
task	qbox		= "ucsclris$src/qbox.cl"

clbye()

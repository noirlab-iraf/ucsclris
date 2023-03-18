#{ Package PREP is Drew Phillips' lris slitmask design software

# cl < "ucsclris$lib/zzsetenv.def"
package	prep, bin = ucsclrisbin$

task	simulator,
	mapmask,
	fabmask,
	gen_igi = "ucsclris$src/x_ucsclris.e"

task	qmask_plot	= "ucsclris$src/qmask_plot.cl"

clbye()

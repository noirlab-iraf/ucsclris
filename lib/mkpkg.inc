# I don't know what this is all about....
# Global MKPKG definitions for the EXAMPLES package.

$set XFLAGS	= "$(XFLAGS) -p ucsclris"
$set XVFLAGS	= "$(XVFLAGS) -p ucsclris"
$set LFLAGS	= "$(LFLAGS) -p ucsclris"

# Uncomment and modify the following to add special file list entries for
# various machine architectures and Fortran compilers.

# $ifeq (MACH, sparc) then
# $include "ucsclris$lib/mkpkg.sf.sun4"
# $else $ifeq (MACH, vms) then
# $include "ucsclris$lib/mkpkg.sf.vms"
# $end

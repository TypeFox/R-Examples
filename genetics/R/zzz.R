## $Id: zzz.R 1330 2007-11-20 20:23:12Z warnes $

# Obsoleted by Proper use the DESCRIPTION and NAMESPACE files
.onAttach <- function(libname, pkgname)
{
	packageStartupMessage("\n")
        packageStartupMessage("NOTE: THIS PACKAGE IS NOW OBSOLETE.\n")
	packageStartupMessage("\n")
	packageStartupMessage("  The R-Genetics project has developed an set of enhanced genetics\n")
	packageStartupMessage("  packages to replace 'genetics'. Please visit the project homepage\n")
        packageStartupMessage("  at http://rgenetics.org for informtion.\n")
	packageStartupMessage("\n")
	
}

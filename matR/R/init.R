
#-----------------------------------------------------------------------------------------
#  package start and finish.
#-----------------------------------------------------------------------------------------

.onAttach <- function (libname, pkgname) {
	packageStartupMessage (tagline ())
	pkgs <- hazPackages()
	if (!all (pkgs)) 
		packageStartupMessage(
			"Suggested package(s) missing: ", 
			collapse (names(pkgs) [!pkgs]))
	}

.Last.lib <- function (libpath) { 
	}

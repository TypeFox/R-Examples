.onAttach <- function(libname, pkgname) {
	packageStartupMessage(paste0("NAEPprimer v", utils::packageDescription("NAEPprimer")$Version, "\n"))
}
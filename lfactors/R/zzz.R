.onAttach <- function(libname, pkgname) {
	packageStartupMessage(paste0("lfactors v", utils::packageDescription("lfactors")$Version, "\n"))
}

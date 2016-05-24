.onAttach <- function(libname, pkgname)
{
ver <- as.character(read.dcf(file.path(libname, pkgname, "DESCRIPTION"), "Version"))
packageStartupMessage("Package epiR ", ver, " is loaded", appendLF = TRUE)
packageStartupMessage("Type help(epi.about) for summary information")
packageStartupMessage("\n")
}
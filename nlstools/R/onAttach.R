.onAttach <- function(libname, pkgname)
{
  packageStartupMessage("\n'nlstools' has been loaded.\n")
  packageStartupMessage("IMPORTANT NOTICE: Most nonlinear regression models and data set examples")
  packageStartupMessage("related to predictive microbiolgy have been moved to the package 'nlsMicrobio'\n")
}

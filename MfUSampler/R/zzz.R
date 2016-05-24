.onAttach <- function(libname, pkgname) {
  RFver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                    fields="Version")
  packageStartupMessage(paste(pkgname, RFver))
  packageStartupMessage("Convenience Functions for Multivariate MCMC Using Univariate Samplers")
  packageStartupMessage("Scientific Computing Group, Sentrana Inc. &\nHeart and Lung Institute, Imperial College London")
}

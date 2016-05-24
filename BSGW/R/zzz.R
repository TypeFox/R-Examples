.onAttach <- function(libname, pkgname) {
  RFver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                    fields="Version")
  packageStartupMessage(paste0("Package: ", pkgname, ", Version: ", RFver))
  packageStartupMessage("Dynamic Survival Model using Bayesian Generalized Weibull Regression.")
  packageStartupMessage("Scientific Computing Group, Sentrana Inc. &\nHeart and Lung Institute, Imperial College London")
}

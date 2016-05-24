.onAttach <- function(libname, pkgname) {
  RFver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                    fields="Version")
  packageStartupMessage(paste0("Package: ", pkgname, ", Version: ", RFver))
  packageStartupMessage("Metropolis-Hastings MCMC using Stochastic Newton Sampler")
  packageStartupMessage("Scientific Computing Group, Sentrana Inc. &")
  packageStartupMessage("Imperial College London")
}

.onAttach <- function(libname, pkgname) {
  RFver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                    fields="Version")
  packageStartupMessage(paste0("Package: ", pkgname, ", Version: ", RFver))
  packageStartupMessage("The R Implementation of 't-walk' MCMC Sampler")
  packageStartupMessage("http://www.cimat.mx/~jac/twalk/")
  packageStartupMessage("For citations, please use:\n")
  packageStartupMessage("Christen JA and Fox C (2010). A general purpose sampling algorithm for")
  packageStartupMessage("continuous distributions (the t-walk). Bayesian Analysis, 5(2),")
  packageStartupMessage("pp. 263-282. <URL:")
  packageStartupMessage("http://ba.stat.cmu.edu/journal/2010/vol05/issue02/christen.pdf>.\n")
}

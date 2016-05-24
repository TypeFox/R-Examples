.onAttach <- function(libname, pkgname) {
  RFver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                    fields="Version")
  packageStartupMessage(paste(pkgname, RFver))
  packageStartupMessage("Stochastic Augmentation of Matched Data Using Restriction Methods")
  packageStartupMessage("Heart and Lung Institute, Imperial College London &\nScientific Computing Group, Sentrana Inc.")
}

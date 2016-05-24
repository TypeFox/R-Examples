.onAttach <- function(libname, pkgname) {
  RFver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                    fields="Version")
  packageStartupMessage(paste0("Package: ", pkgname, ", Version: ", RFver))
  packageStartupMessage("Cause-specific Framework for Competing-risk survival analysis")
  packageStartupMessage("School of Public Heath, Imperial College London &")
  packageStartupMessage("Scientific Computing Group, Sentrana Inc.")
}

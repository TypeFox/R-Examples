.onAttach <- function(libname, pkgname){
  packageStartupMessage(paste("\nmpm version ", packageDescription("mpm")$Version, 
          "\n", sep = ""))
}


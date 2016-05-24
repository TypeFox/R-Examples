#--------------------------------------------------------------------
.onAttach <- function(libname, pkgname){
#--------------------------------------------------------------------
#   pkg.info <- utils::packageDescription(pkgname, libname, fields = c("Title", "Version", "Date"))
    pkg.info <- drop( read.dcf( file = system.file("DESCRIPTION", package = "qcr"),
                      fields = c("Title", "Version", "Date") ))
    packageStartupMessage( 
      paste("\n Package qcr:", pkg.info["Title"], "\n"),
      paste(" version ", pkg.info["Version"], " (built on ", pkg.info["Date"], ").\n", sep=""),
      " Copyright Miguel A. Flores Sanchez 2014. \n")
}

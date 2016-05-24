#--------------------------------------------------------------------
.onAttach <- function(libname, pkgname){
#--------------------------------------------------------------------
#   pkg.info <- utils::packageDescription(pkgname, libname, fields = c("Title", "Version", "Date"))
    pkg.info <- drop( read.dcf( file = system.file("DESCRIPTION", package = "mme"),
                      fields = c("Title", "Version", "Date") ))
    packageStartupMessage( 
      paste("\n Package mme:", pkg.info["Title"], "\n"),
      paste(" Version ", pkg.info["Version"], " (built on ", pkg.info["Date"], ") is now loaded.\n", sep=""),
      " Copyright E. Lopez-Vizcaino, M.J. Lombardia and D. Morales \n")
}

#--------------------------------------------------------------------
#   zzz.R (npsp package)
#--------------------------------------------------------------------
#   .onAttach(libname, pkgname)
#
#   (c) R. Fernandez-Casal         Last revision: Apr 2013
#--------------------------------------------------------------------


#--------------------------------------------------------------------
.onAttach <- function(libname, pkgname){
#--------------------------------------------------------------------
#   pkg.info <- utils::packageDescription(pkgname, libname, fields = c("Title", "Version", "Date"))
    pkg.info <- drop( read.dcf( file = system.file("DESCRIPTION", package = "npsp"),
                      fields = c("Title", "Version", "Date") ))
    packageStartupMessage( 
      paste(" Package npsp:", pkg.info["Title"], "\n"),
      paste(" version ", pkg.info["Version"], " (built on ", pkg.info["Date"], ").\n", sep=""),
      " Copyright R. Fernandez-Casal 2012-2014.\n",      
      " Type demo(package = 'npsp') to obtain the list of available demos.\n")
}

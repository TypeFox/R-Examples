is.fdPar <- function(fdParobj) {
#  check whether FDPAROBJ is a functional data object

#  Last modified 20 November 2005

  if (inherits(fdParobj, "fdPar")) return(TRUE) else return(FALSE)
}

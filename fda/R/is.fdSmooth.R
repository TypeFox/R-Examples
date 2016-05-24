is.fdSmooth <- function(fdSmoothobj) {
#  check whether FDPAROBJ is a functional data object

#  Last modified 20 November 2005

  if (inherits(fdSmoothobj, "fdSmooth")) return(TRUE) else return(FALSE)
}
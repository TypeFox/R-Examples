is.fd <- function(fdobj) {
#  check whether FDOBJ is a functional data object

#  Last modified 20 November 2005

  if (inherits(fdobj, "fd")) return(TRUE) else return(FALSE)
}

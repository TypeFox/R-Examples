is.basis <- function(basisobj) {
#  check whether BASISOBJ is a functional data basis object

#  Last modified 20 November 2005

  if (inherits(basisobj, "basisfd")) return(TRUE) else return(FALSE)
}

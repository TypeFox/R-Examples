getbasisrange <- function(basisobj){
#  GETBASISRANGE   Extracts the range from basis object BASISOBJ.
#  R port 2007.11.28 by Spencer Graves  
#  previously modified 30 June 1998

  if(!is.basis(basisobj))
    stop("'basisobj' is not a functional basis object.")
#
  rangeval <- basisobj$rangeval
  rangeval
}

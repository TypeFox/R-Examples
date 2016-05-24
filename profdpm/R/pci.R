pci <- function( x1, x2 ) {
  ###################################################
  #do some argument checking
  if( length( x1) != length( x2 ) ) { stop("x1 and x2 lengths differ") }
  ret <- .Call("pci", as.factor(x2), as.factor(x1), PACKAGE="profdpm")
  return( ret )
}

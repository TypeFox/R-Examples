# Cleaning a matrix from nas by row
# If X is not a matrix or a vector it will not be cleaned
# If there are no nas we return list list("cleanX" = X, "banned" = numeric(), "nBanned" = 0)
# If there are nas we return the same list, where x has been cleaned

.clean <- function(X, verbose = FALSE){
  
  if( is.vector(X, mode = "numeric") ) X <- as.matrix(X, length(X), 1)
  
  if( is.matrix(X) )
  {
    out <- .Call( "cleanStats", inMat = X, PACKAGE = "synlik" )
    
    if( out$nBanned > 0)
    {
      if( verbose && (out$nBanned > 0.75 * nrow(X)) ) 
        warning( paste(out$nBanned, "out of", nrow(X), "statistics vectors", "contain NAs and will not be used") )
      
      return( out )
    }
    
  } 
  
  return( list("nBanned" = 0, "banned" = numeric(), "cleanX" = X) )
  
}

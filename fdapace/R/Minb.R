# This function is used to find the minimum bandwidth choice
# where local window contains at least "numPoints" points
# Input x  : n x 1 vector
# Input numPoints: an integer specifying the number of points in a local window
# for local weighted constant, numPoints is at least 1
# for local weighted linear, numPoints is at least 2
# Output b: the minimum bandwidth choice for vector x

Minb <- function(x, numPoints){ 
  x = sort(unique(x));     # Unique is added to ensure that we do not have a degenerate design
  n = length(x);
  if( (numPoints<1) || (numPoints > n) ){
   warning("Invalid number of minimum points specified\n")
   return(NaN)
  }
  if(numPoints > 1){ 
    return( max(x[numPoints:n]-x[1:(n-numPoints+1)]) )
  }else{ 
    return( max(   (x[2:n]-x[1:(n-1)])/2)   )
  } 
}

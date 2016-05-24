#' Spherical Coordinates Rotation
#'
#' Rotate to the direction of interest in spherical coordinates by contrasts (e.g. 1, -1, -1).
#' @usage rotate3d(T, contrast)
#' @param T a numeric matrix (normally z-score converted) with 3 columns. The rows are genes or phosphorylation
#' sites and the columns are treatments vs control statistics.
#' @param contrast the direction to be tested for enrichment. Specified as contrast in a triplet (see example).
#' @return A rotated matrix with respect to the direction of interest.
#' @export 
#' @examples
#' # load the example data
#' data(PM)
#' 
#' # convert statistics into z-scores
#' PM.zscores <- apply(PM, 2, function(x){qnorm(rank(x)/(nrow(PM)+1))})
#' 
#' # Rotate the matrix by contrast 1, -1, -1 (i.e. up-regulation, down-regulation, dow-regulation).
#' PM.rotated <- rotate3d(PM.zscores, contrast = c(1, -1, -1))
#' 
rotate3d = function(T, contrast = c(1,1,1)){

  a = contrast
  b= c(1,1,1) # the original direction in contrast form is c(1,1,1)
  
  if(sum(a == b)==3)return(T)
  if(sum(a == -b)==3)return(-T)
  
  a = a/sqrt(sum(a^2))
  b = b/sqrt(sum(b^2))
  u = c((a[2]*b[3] - a[3]*b[2]) , (a[3]*b[1] - a[1]*b[3]) ,(a[1]*b[2] - a[2]*b[1]) )
  u = u/sqrt(sum(u^2))
  R = matrix(0,3,3)
  the = acos(sum(a*b))
  R[1,1] = cos(the) + u[1]^2*(1-cos(the))
  R[2,2] = cos(the) + u[2]^2*(1-cos(the))
  R[3,3] = cos(the) + u[3]^2*(1-cos(the))           
  R[1,2] = u[1]*u[2]*(1-cos(the)) - u[3]*sin(the)                 
  R[1,3] = u[1]*u[3]*(1-cos(the)) + u[2]*sin(the)                 
  R[2,1] = u[2]*u[1]*(1-cos(the)) + u[3]*sin(the)                 
  R[2,3] = u[2]*u[3]*(1-cos(the)) - u[1]*sin(the)                 
  R[3,1] = u[3]*u[1]*(1-cos(the)) - u[2]*sin(the)                 
  R[3,2] = u[3]*u[2]*(1-cos(the)) + u[1]*sin(the)                 
  return(t(R%*%t(T)))
}

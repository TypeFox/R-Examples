#' Polar Coordinates Rotation
#'
#' Rotate to the direction of interest in polar coordinates by degree (e.g. pi/4).
#' @usage rotate2d(T, degree = 0)
#' @param T a numeric matrix (normally z-score converted) with 2 columns. The rows are genes or phosphorylation
#' sites and the columns are treatments vs control statistics.
#' @param degree the direction to be tested for enrichment. Specified as degree from to the original direction.
#' @return A rotated matrix with respect to the direction of interest.
#' @export 
#' @examples
#' # load the phosphoproteomics dataset
#' data(HEK)
#' 
#' # convert statistics into z-scores
#' HEK.zscores <- apply(HEK, 2, function(x){qnorm(rank(x)/(nrow(HEK)+1))})
#' 
#' # Rotate the matrix by 1/2 pi (i.e. down-regulation, dow-regulation).
#' HEK.rotated <- rotate2d(HEK.zscores, degree = pi/2)
#' 
rotate2d = function(T, degree = 0){
  theta = atan2(T[,2], T[,1]) + degree  
  r = sqrt(T[,1]^2 + T[,2]^2)
  x2 = r*cos(theta)
  y2 = r*sin(theta)
  return(cbind(x2, y2))
}

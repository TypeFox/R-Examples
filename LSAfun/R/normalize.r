#' @export normalize

normalize <- function(x){
  
  norm <- sqrt(sum(x^2))
  
  normvec <- x/norm
  normvec
}
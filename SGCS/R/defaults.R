#' Default range vector
#' @param x point pattern, internal format
#' @param r most likely missing. if not , returned.
#' Compute the range vector
#'  

default_r <- function(x, r){
  if(!missing(r)) return(r)
  r <- min(apply(x$bbox, 2, diff))/3
  rlarge <- sqrt(1000/(pi * x$n/x$area))
  m <- min(rlarge, r)
  seq(0, m, length=50)
}
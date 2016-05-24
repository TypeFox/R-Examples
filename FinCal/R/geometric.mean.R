#' Geometric mean return
#' 
#' @param r returns over multiple periods
#' @export
#' @examples
#' geometric.mean(r=c(-0.0934, 0.2345, 0.0892))
geometric.mean <- function(r){
  rs <- r + 1
  return(prod(rs)^(1/length(rs))-1)
}

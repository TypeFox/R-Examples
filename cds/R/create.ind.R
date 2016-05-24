#' Create Indicator Matrix
#' 
#' Create an indicator matrix.
#' 
#' @param grp A grouping vector.
#' @keywords utility
#' @export create.ind
create.ind <- function(grp) {
  K <- length(unique(grp))
  n <- length(grp)
  mat <- matrix(0, nrow = n, ncol = K)
  out <- apply(rbind(1:K, mat), 2, function(x, a) x[-1] + (a == x[1]), a = grp)
  out
}

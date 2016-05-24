#' @rdname runif_sphere
#' @export
#' @examples
#' 
#' set.seed(1234)
#' runif_sphere(n = 10)
runif_cube <- function (n, m = 2, max.dist = 1, at = FALSE, nr.dist = 21) { 
  if(at){
    spec.dist <- seq(from = 0, to = max.dist, length = nr.dist)
#     nper <- ceiling(n/nr.dist)  
    nper <- ceiling(n * spec.dist / sum(spec.dist))
    nper[1] <- 1
    n <- sum(nper)
    out <- matrix(runif(n*m, min = -max.dist, max = max.dist), nrow = n, 
                  ncol = m)
    spec.dist <- rep(spec.dist, nper)
    out <- sweep(out, MARGIN = 1, STATS = apply(out, 1, function(x) max(abs(x))),
                 FUN = "/")
    out <- spec.dist*out
  } else out <- matrix(runif(n*m, min = -max.dist, max = max.dist), nrow = n, 
                       ncol = m)
  out
}

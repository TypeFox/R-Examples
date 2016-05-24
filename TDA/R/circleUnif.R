circleUnif <-
function(n, r = 1){

  if (!is.numeric(n) || length(n) != 1 || n < 0) {
    stop("n should be a nonnegative integer")
  }
  if (!is.numeric(r) || length(r) != 1 || r < 0) {
    stop("r should be a nonnegative number")
  }

  ########################
  th <- stats::runif(n, 0, 2 * pi)
  x1 <- r * cos(th)
  x2 <- r * sin(th)
  X <- cbind(x1, x2)
  return(X)
}

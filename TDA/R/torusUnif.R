torusUnif <-
function(n, a, c) {

  if (!is.numeric(n) || length(n) != 1 || n < 0) {
    stop("n should be a nonnegative integer")
  }
  if (!is.numeric(a) || length(a) != 1 || a < 0) {
    stop("a should be a nonnegative number")
  }
  if (!is.numeric(c) || length(c) != 1 || c < 0) {
    stop("c should be a nonnegative number")
  }

  n <- floor(n)
  theta <- NULL
  while (length(theta) < n){
    xvec <- stats::runif(1, 0, 2 * pi)
    yvec <- stats::runif(1, 0, 1 / pi)
    fx <- (1 + (a / c) * cos(xvec)) / (2 * pi)
    if (yvec < fx) {
      theta <- c(theta, xvec)
    }
  }

  phi <- stats::runif(n, 0, 2 * pi)
  x <- (c + a * cos(theta)) * cos(phi)
  y <- (c + a * cos(theta)) * sin(phi)
  z <- a * sin(theta)
  
  out <- cbind(x, y, z)
  return(out)
}

sphereUnif <-
function(n, d, r = 1) {

  if (!is.numeric(n) || length(n) != 1 || n < 0) {
    stop("n should be a nonnegative integer")
  }
  if (!is.numeric(d) || length(d) != 1 || d < 0) {
    stop("d should be a nonnegative integer")
  }
  if (!is.numeric(r) || length(r) != 1 || r < 0) {
    stop("r should be a nonnegative number")
  }
     
  ########################
  X <- array(stats::rnorm(n * (d+1)), dim = c(n, d + 1))
  for (i in 1:n) {
    while(norm(X[i, ], "2") == 0) {
      X[i, ] <- stats::rnorm(d + 1)
    }
    X[i,] <- r * X[i, ] / norm(X[i, ], "2")
  }
  return(X)     
}

rmLaplace <- 
function(n = 1, center = rep(0, nrow(Scatter)), Scatter = diag(length(center)))
{
  if (n <= 0)
    stop("n must be a positive integer")
  if (length(center) != nrow(Scatter))
    stop("center and scatter have non-conforming size")
  p <- nrow(Scatter)

  y <- matrix(0, nrow = n, ncol = p)
  dy <- dim(y)
  # call C code
  y <- .C("rand_laplace", 
          y = as.double(y),
          dims = as.integer(dy),
          center = as.double(center),
          Scatter = as.double(Scatter))$y
  y <- matrix(y, nrow = n, byrow = TRUE)
  y
}

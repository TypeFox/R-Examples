rmt <- 
function(n = 1, mean = rep(0, nrow(Sigma)), Sigma = diag(length(mean)), eta = .25)
{
  if (n <= 0)
    stop("n must be a positive integer")
  if (length(mean) != nrow(Sigma))
    stop("mean and Sigma have non-conforming sizes")
  if ((eta < 0 ) || (eta >= 1/2))
      stop("eta must be in [0,1/2)")
  p <- nrow(Sigma)

  y <- matrix(0, nrow = n, ncol = p)
  dy <- dim(y)
  # call C code
  if (eta == 0) {
    y <- .C("rand_norm",
            y = as.double(y),
            dims = as.integer(dy),
            mean = as.double(mean),
            Sigma = as.double(Sigma))$y
  }
  else {
    y <- .C("rand_student",
            y = as.double(y),
            dims = as.integer(dy),
            mean = as.double(mean),
            Sigma = as.double(Sigma),
            eta = as.double(eta))$y
  }
  y <- matrix(y, nrow = n, byrow = TRUE)
  y
}

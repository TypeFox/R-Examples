rmvnorm <- function (n = 1, mean, var) {
  p <- length(mean)
  if (!all(dim(var)==c(p,p)))
    stop("length of ",sQuote("mean")," must equal the dimension of the square matrix ",sQuote("var"))
  cf <- t(chol(var))
  matrix(mean,p,n)+cf%*%matrix(rnorm(p*n),p,n)
}


################################
#### Spatial median regression
#### Tsagris Michail 10/2014
#### Biman Chakraborty (2003) On multivariate quantile regression
#### Journal of Statistical Planning and Inference
#### http://www.stat.nus.edu.sg/export/sites/dsap/research/documents/tr01_2000.pdf
#### mtsagris@yahoo.gr
################################

spatmed.reg <- function(y, x, xnew = NULL) {
  ## y contains the dependent variables
  ## x contains the independent variable(s)

  runtime <- proc.time()
  y <- as.matrix(y)
  x <- as.matrix(x)
  d <- ncol(y)  ## dimensionality of y
  x <- cbind(1, x)  ## add the constant term
  p <- ncol(x)  ## dimensionality of x
  z <- list(y = y, x = x)
  ## medi is the function to perform median regression
  medi <- function(beta, z) {
    y <- z$y
    x <- z$x
    p <- ncol(x)
    be <- matrix(beta, nrow = p)
    est <- x %*% be
    sum( sqrt( rowSums((y - est)^2) ) )
  }
  ## we use nlm and optim to obtain the beta coefficients
  ini <- matrix(nrow = p, ncol = d)
  for (i in 1:d)  ini[, i] <- coef( quantreg::rq(y[, i] ~ x[, -1]) )
  ini <- as.vector(ini)
  qa <- nlm(medi, ini, z = z, iterlim = 10000)
  qa <- optim(qa$estimate, medi, z = z, control = list(maxit = 20000))
  qa <- optim(qa$par, medi, z = z, control = list(maxit = 20000),
  hessian = TRUE)
  beta <- matrix( qa$par, ncol = ncol(y) )

  if ( is.null(xnew) ) {
    est = x %*% beta
  } else {
    xnew <- cbind(1, xnew)
    xnew <- as.matrix(xnew)
    est <- xnew %*% beta
  }

  seb <- sqrt( diag( solve(qa$hessian) ) )
  seb <- matrix( seb, ncol = ncol(y) )

  if ( is.null(colnames(y)) ) {
    colnames(seb) <- colnames(beta) <- paste("Y", 1:d, sep = "")
  } else  colnames(seb) <- colnames(beta) <- colnames(y)

  if ( is.null(colnames(x)) ) {
    p <- ncol(x) - 1
    rownames(beta) <- c("constant", paste("X", 1:p, sep = "") )
    if ( !is.null(seb) )  rownames(seb) <- c("constant", paste("X", 1:p, sep = "") )
  } else {
    rownames(beta)  <- c("constant", colnames(x)[-1] )
    if  ( !is.null(seb) ) rownames(seb) <- c("constant", colnames(x)[-1] )
  }

  if (d == 1)  est <- as.vector(est)
  runtime <- proc.time() - runtime

  list(runtime = runtime, beta = beta, seb = seb, est = est)
}

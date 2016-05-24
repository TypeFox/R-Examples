################################
#### Multivariate (and univariate) ridge regression
#### Tsagris Michail 10/2014
#### mtsagris@yahoo.gr
################################

ridge.reg <- function(y, x, lambda, B = 1, xnew = NULL) {
  ## xnew is the new independent variables values
  ## whose values of y you want to estimate
  ## by default xnew is the x, so you will get the fitted values
  ## y is a real valued vector
  ## x contains the independent variable(s)
  ## lambda is the ridge regularization parameter
  ## if lambda=0, the classical multivariate regression is implemented
  ## B is for bootstrap estimation of the standard errors of the betas
  ## if pred is TRUE it means that you want to predict new y
  ## but if xnew is x (by default), the pred is not important
  ## the pred is important if xnew is not x
  y <- as.vector(y)
  if ( all( y > 0 & y < 1 ) ){
    y <- log(y / ( 1 - y) ) ## logistic normal
  }
  x <- as.matrix(x)
  n <- length(y)  ## sample size
  p <- ncol(x)  ## dimensionality of x
  my <- mean(y)
  yy <- y - my ## center the dependent variables
  xx <- scale(x)[1:n, ]  ## standardize the independent variables
  W <- solve( crossprod(xx) + lambda * diag(p) )
  beta <- W %*% crossprod(xx, yy)
  est <- as.vector( xx %*% beta + my )
  va <- var(y - est) * (n - 1) / (n - p - 1)
  vab <- kronecker(va, W %*% crossprod(xx) %*% W  )
  seb <- matrix( sqrt( diag(vab) ), nrow = p )
  if (B > 1) { ## bootstrap estimation of the standard errors
    be <- matrix(nrow = B, ncol = p )
    for ( i in 1:B) {
      id <- sample(1:n, n, replace = TRUE)
      yb <- yy[id, ]  ;  xb <- xx[id, ]
      W <- solve( crossprod(xb) + lambda * diag(p) )
      be[i, ] <- W %*% crossprod(xb, yb)
    }
    seb <- matrix( apply(be, 2, sd), nrow = p ) ## bootstrap standard errors of betas
  }
  ## seb contains the standard errors of the coefficients
  if ( is.null( colnames(x) ) ) {
    names(seb) <- paste("X", 1:p, sep = "")
    names(beta) <- paste("X", 1:p, sep = "")
  } else  names(seb) <- names(beta) <- colnames(x)
  if ( !is.null(xnew) ) {
    mx <- colMeans(x)
    s <- apply(x, 2, sd)
    s <- diag(1/s)
    xnew <- as.matrix(xnew)
    xnew <- matrix(xnew, ncol = p)
    nu <- nrow(xnew)
    xnew <- ( xnew - rep( mx, rep(nu, p) ) ) %*% s ## scale the xnew values
    est <- xnew %*% beta + my
  } else est <- x %*% beta + my
  est <- as.vector(est)
  list(beta = beta, seb = seb, est = est)
}

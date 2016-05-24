################################
#### alfa-regression
#### Tsagris Michail 11/2015
#### mtsagris@yahoo.gr
#### References: Tsagris Michail (2013)
#### Regression analysis with compositional data containing zero values
#### Chilean Journal of Statistics, 6(2): 47-57
################################

alfa.reg <- function(y, x, a, xnew = NULL, yb = NULL) {
  ## y is the compositional data (dependent variable)
  ## x is the independent variables
  ## a is the value of alpha
  y <- as.matrix(y)
  y <- y/rowSums(y)  ## makes sure y is compositional data
  x <- as.matrix(x)
  p <- ncol(x)   ;   n <- nrow(x)

  if ( !is.null(xnew) ) {
    ## if the xnew is the same as the x, the classical fitted values
    ## will be returned. Otherwise, the estimated values for the
    ## new x values will be returned.
    if (p == 1) {
      mx <- mean(x)
      s <- sd(x)
      xnew <- (xnew - mx) / s
    } else {
      x <- as.matrix(x)
      mx <- colMeans(x)
      s <- apply(x, 2, sd)
      s <- diag(1/s)
      xnew <- as.matrix(xnew)
      xnew <- matrix(xnew, ncol = p)
      nu <- nrow(xnew)
      xnew <- ( xnew - rep(mx, rep(nu, p)) ) %% s  ## standardize the xnew values
    }
    xnew <- cbind(1, xnew)
  }

  x <- scale(x)[1:n, ]  ## standardize the independent variables
  x <- as.matrix( cbind(1, x) )
  d <- ncol(y) - 1  ## dimensionality of the simplex

  ## internal function for the alfa-regression
  reg <- function(para, z){
    ya <- z$ya   ;   x <- z$x
    d <- ncol(ya)  ;  p <- ncol(x)
    m0 <- numeric(d)
    n <- nrow(ya)
    be <- matrix(para, byrow = TRUE, ncol = d)
    mu1 <- cbind( 1, exp(x %*% be) )
    mu <- mu1 / rowSums(mu1)
    ma <- alfa(mu, a)$aff
    esa <- ya - ma
    sa <- crossprod(esa) / (n - p)
    f <- ( n/2 ) * log( det(sa) ) + 0.5 * sum( mahalanobis(esa, m0, sa) )
    f
  }

  if ( a == 0 ) {
    mod <- comp.reg(y, x[, -1], yb = yb)
    beta <- mod$beta
    seb <- mod$seb
    runtime <- mod$runtime

  } else {
    runtime <- proc.time()

    if ( is.null(yb) ) {
      yb <- alfa(y, a)$aff
    } else {
      yb <- yb
    }

    z <- list(ya = yb, x = x)
    ini <- as.vector( coef(lm(yb ~ x[, -1])) )
    qa <- nlminb( ini, reg, z = z, control = list(iter.max = 1000) )
    qa <- optim( qa$par, reg, z = z, control = list(maxit = 5000) )
    qa <- optim( qa$par, reg, z = z, control = list(maxit = 5000) )
    qa <- optim( qa$par, reg, z = z, control = list(maxit = 5000) )
    qa <- optim( qa$par, reg, z = z, control = list(maxit = 5000), hessian = TRUE )
    beta <- matrix(qa$par, byrow = TRUE, ncol = d)
    seb <- sqrt( diag( solve( qa$hessian) ) )
    seb <- matrix(seb, byrow = TRUE, ncol = d)

    runtime <- proc.time() - runtime
  }

  if ( !is.null(xnew) ) {
    mu <- cbind( 1, exp(xnew %*% beta) )
    est <- mu/rowSums(mu)
  } else {
    mu <- cbind(1, exp(x %*% beta) )
    est <- mu/rowSums(mu)
  }

  if ( is.null(colnames(x)) ) {
    p <- ncol(x) - 1
    rownames(beta) <- c("constant", paste("X", 1:p, sep = "") )
    if ( !is.null(seb) )  rownames(seb) <- c("constant", paste("X", 1:p, sep = "") )
  } else {
    rownames(beta)  <- c("constant", colnames(x)[-1] )
    if  ( !is.null(seb) ) rownames(seb) <- c("constant", colnames(x)[-1] )
  }

  list(runtime = runtime, beta = beta, seb = seb, est = est)
}

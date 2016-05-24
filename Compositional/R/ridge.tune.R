################################
#### Ridge regression tuning of the
#### lambda parameter via k-fold cross validation
#### Tsagris Michail 1/2016
#### mtsagris@yahoo.gr
################################

ridge.tune <- function(y, x, M = 10, lambda = seq(0, 2, by = 0.1),
  mat = NULL, ncores = 1, graph = FALSE) {

  ## y is the univariate or multivariate dependent variable
  ## x contains the independent variables(s)
  ## M is the number of folds, set to 10 by default
  ## lambda is a vector with a grid of values of lambda
  ## ncores is the number of cores to use

  y <- as.matrix(y)  ## makes sure y is a matrix
  x <- as.matrix(x)
  n <- nrow(y)  ## sample size
  k <- length(lambda)

  di <- ncol(y)  ## dimensionality of y
  p <- ncol(x)  ## dimensionality of x
  if ( is.null(mat) ) {
    nu <- sample(1:n, min( n, round(n / M) * M ) )
    ## It may be the case this new nu is not exactly the same
    ## as the one specified by the user
    ## to a matrix a warning message should appear
    options(warn = -1)
    mat <- matrix( nu, ncol = M ) # if the length of nu does not fit
  } else  mat <- mat

  M <- ncol(mat)
  rmat <- nrow(mat)
  msp <- matrix( nrow = M, ncol = k)

  ## deigma will contain the positions of the test set
  ## this is stored but not showed in the end
  ## the user can access it though by running
  ## the commands outside this function

  if (ncores == 1) {
    runtime <- proc.time()
    for (vim in 1:M) {
      ytest <- as.matrix( y[mat[, vim], ] )  ## test set dependent vars
      ytrain <- as.matrix( y[-mat[, vim], ] )  ## train set dependent vars
      my <- colMeans(ytrain)
      yy <- ytrain - rep(my, rep(rmat, di) )  ## center the dependent variables
      xtrain <- as.matrix( x[ -mat[, vim], ] )  ## train set independent vars
      mx <- colMeans(xtrain)
      xtest <- as.matrix( x[ mat[, vim], ] )  ## test set independent vars
      s <- apply(xtrain, 2, sd)
      s <- diag(1/s)
      xtest <- ( xtest - rep( mx, rep(rmat, p) ) ) %*% s ## standardize the xtest
      xx <- scale(xtrain)[1:c(n - rmat), ]  ## standardize the independent variables
      sa <- svd(xx)
      u <- sa$u   ;   d <- sa$d   ;   v <- sa$v
      for ( i in 1:k ) {
        beta <- ( v %*% diag( d / ( d^2 + lambda[i] ) ) %*% t(u) ) %*% yy
        est <- xtest %*% beta + my
        msp[vim, i] <- sum( (ytest - est)^2 ) / rmat
      }
    }
    runtime <- proc.time() - runtime

  } else {
    runtime <- proc.time()
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    pe <- numeric(k)
    msp <- foreach(vim = 1:M, .combine = rbind) %dopar% {
      ytest <- as.matrix( y[mat[, vim], ] )  ## test set dependent vars
      ytrain <- as.matrix( y[-mat[, vim], ] )  ## train set dependent vars
      my <- colMeans(ytrain)
      yy <- ytrain - rep(my, rep(rmat, di) )  ## center the dependent variables
      xtrain <- as.matrix( x[ -mat[, vim], ] )  ## train set independent vars
      mx <- colMeans(xtrain)
      xtest <- as.matrix( x[ mat[, vim], ] )  ## test set independent vars
      s <- apply(xtrain, 2, sd)
      s <- diag(1/s)
      xtest <- ( xtest - rep( mx, rep(rmat, p) ) ) %*% s ## standardize the xtest
      xx <- scale(xtrain)[1:c(n - rmat), ]  ## standardize the independent variables
      sa <- svd(xx)
      u <- sa$u   ;   d <- sa$d   ;   v <- sa$v
      for ( i in 1:k ) {
        beta <- ( v %*% diag( d / ( d^2 + lambda[i] ) ) %*% t(u) ) %*% yy
        est <- xtest %*% beta + my
        pe[i] <- sum( (ytest - est)^2 ) / rmat
      }
      return(pe)
    }
    runtime <- proc.time() - runtime
    stopCluster(cl)
  }
  mspe <- colMeans(msp)
  bias <- msp[ , which.min(mspe)] - apply(msp, 1, min)  ## TT estimate of bias
  estb <- mean( bias )  ## TT estimate of bias

  if (graph == TRUE) {
    plot(lambda, mspe, type = 'b', ylim = c(min(mspe), max(mspe)),
    ylab = "Mean squared error of prediction",
    xlab = expression(paste(lambda, " values")) )
  }

  names(mspe) <- lambda
  performance <- c( min(mspe) + estb, estb)
  names(performance) <- c("MSPE", "Estimated bias")
  list(msp = msp, mspe = mspe, lambda = which.min(mspe), performance = performance,
  runtime = runtime)
}

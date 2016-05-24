################################
#### Principal components regression for binary and poisson regression
#### Selection of the number of principal components
#### via K-fold cross validation
#### Tsagris Michail 1/2016
#### mtsagris@yahoo.gr
################################

glmpcr.tune <- function(y, x, M = 10, maxk = 10, mat = NULL,
                        ncores = 1, graph = TRUE) {
  ## y is the UNIVARIATE dependent variable
  ## y is either a binary variable (binary logistic regression)
  ## or a discrete variable (Poisson regression)
  ## x contains the independent variables
  ## fraction denotes the percentage of observations
  ## to be used as the test set
  ## the 1-fraction proportion of the data will be the training set
  ## R is the number of cross validations
  ## if ncores==1, then 1 processor is used, otherwise more are
  ## used (parallel computing)
  x <- as.matrix(x)
  y <- as.vector(y)
  n <- nrow(x)
  p <- ncol(x)
  if ( maxk > p ) maxk <- p  ## just a check

  if ( is.null(mat) ) {
    nu <- sample(1:n, min( n, round(n / M) * M ) )
    ## It may be the case this new nu is not exactly the same
    ## as the one specified by the user
    ## to a matrix a warning message should appear
    options(warn = -1)
    mat <- matrix( nu, ncol = M )
  } else  mat <- mat

  M <- ncol(mat)
  rmat <- nrow(mat)
  ntrain = n - rmat
  msp <- matrix( nrow = M, ncol = maxk )

  ## deigma will contain the positions of the test set
  ## this is stored but not showed in the end
  ## the user can access it though by running
  ## the commands outside this function
  rmat <- nrow(mat)
  ntrain <- n - rmat

  if ( length(unique(y)) == 2 ) {
    oiko <- "binomial"
  } else oiko <- "poisson"

  if (ncores == 1) {
    runtime <- proc.time()
    for (vim in 1:M) {
      ytest <- as.vector( y[mat[, vim] ] )  ## test set dependent vars
      ytrain <- as.vector( y[-mat[, vim] ] )  ## train set dependent vars
      xtrain <- as.matrix( x[-mat[, vim], ] )  ## train set independent vars
      xtest <- as.matrix( x[mat[, vim], ] )  ## test set independent vars

      mx <- colMeans(xtrain)
      s <- apply(xtrain, 2, sd)
      s <- diag(1/s)
      xtrain <- scale(xtrain)[1:ntrain, ]  ## standardize the independent variables
      eig <- eigen( crossprod(xtrain) )  ## eigen analysis of the design matrix
      vec <- eig$vectors  ## eigenvectors, or principal components
      z <- xtrain %*% vec  ## PCA scores
      xnew <- ( xtest - rep(mx, rep(rmat, p)) ) %*% s ## standardize the xnew values

      for ( j in 1:maxk) {
        mod <- glm(ytrain ~ z[, 1:j], family = oiko )
        b <- coef(mod)
        be <- vec[, 1:j] %*% as.matrix( b[-1] )
        es <- as.vector( xnew %*% be ) + b[1]

        if (oiko == "binomial") {
          est <- as.vector(  exp(es)/(1 + exp(es))  )
          ri <-  -2 *( ytest * log(est) + (1 - ytest) * log(1 - est) )
        } else {
          est <- as.vector( exp(es) )
          ri <- 2 * ( ytest * log(ytest / est) )
        }
        msp[vim, j] <- sum( ri )
      }
    }
    runtime <- proc.time() - runtime

  } else {
    runtime <- proc.time()
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    er <- numeric(maxk)
    msp <- foreach(vim = 1:M, .combine = rbind) %dopar% {
      ## will always be the same
      ytest <- as.vector( y[mat[, vim] ] )  ## test set dependent vars
      ytrain <- as.vector( y[-mat[, vim] ] )  ## train set dependent vars
      xtrain <- as.matrix( x[-mat[, vim], ] )  ## train set independent vars
      xtest <- as.matrix( x[mat[, vim], ] )  ## test set independent vars

      mx <- colMeans(xtrain)
      s <- apply(xtrain, 2, sd)
      s <- diag(1/s)
      xtrain <- scale(xtrain)[1:ntrain, ]  ## standardize the independent variables
      eig <- eigen( crossprod(xtrain) )  ## eigen analysis of the design matrix
      vec <- eig$vectors  ## eigenvectors, or principal components
      z <- xtrain %*% vec  ## PCA scores
      xnew <- ( xtest - rep(mx, rep(rmat, p)) ) %*% s ## standardize the xnew values
      xnew <- cbind(1, xnew %*% vec)

      for ( j in 1:maxk) {
        mod <- glm(ytrain ~ z[, 1:j], family = oiko )
        b <- coef(mod)
        be <- vec[, 1:j] %*% as.matrix( b[-1] )
        es <- as.vector( xnew %*% be ) + b[1]

        if (oiko == "binomial") {
          est <- as.vector(  exp(es)/(1 + exp(es))  )
          ri <-  -2 *( ytest * log(est) + (1 - ytest) * log(1 - est) )
        } else {
          est <- as.vector( exp(es) )
          ri <- 2 * ( ytest * log(ytest / est) )
        }
        er[j] <- sum( ri )
      }
      return(er)
    }
    stopCluster(cl)
    runtime <- proc.time() - runtime
  }

  mpd <- colMeans(msp)
  bias <- msp[ ,which.min(mpd)] - apply(msp, 1, min)  ## TT estimate of bias
  estb <- mean( bias )  ## TT estimate of bias

  if (graph == TRUE) {
    plot(1:maxk, mpd, xlab = "Number of principal components",
         ylab = "Mean predicted deviance", type = "b")
  }

  names(mpd) <- paste("PC", 1:maxk, sep = " ")
  performance <- c( min(mpd) + estb, estb)
  names(performance) <- c("MPD", "Estimated bias")
  list(msp = msp, mpd = mpd, k = which.min(mpd), performance = performance, runtime = runtime)
}

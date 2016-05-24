knnreg.tune <- function(y, x, M = 10, A = 10, ncores = 1, res = "eucl",
  type = "euclidean", estim = "arithmetic", mat = NULL, graph = FALSE) {

  ## y is the multivariate (or univariate) dependent variable
  ## x contains the independent variables(s)
  ## M is the number of folds, set to 10 by default
  ## it is assumed that the training set contains at least 11 observations
  ## A is the highest number of nearest neighbours
  ## res is the type of response, Euclidean ("eucl") or spherical ("spher")
  ## ncores specifies how many cores to use
  ## type is for the distance, Euclidean or Manhattan distance.
  ## The function dist()  allows for more distance types
  ## which can of course use here.
  ## Type ?dist so see more
  ## estim is either 'arithmetic', 'harmonic'. How to calculate the
  ## estimated value of the Y using the arithmetic mean or the
  ## harmonic mean of the closest observations.

  y <- as.matrix(y)
  x <- as.matrix(x)
  n <- nrow(y)
  d <- ncol(y)
  ina <- 1:n

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
  per <- matrix(nrow = M, ncol = A - 1)

  if (type == "spher") {

    runtime <- proc.time()
    x <- x / sqrt( rowSums(x^2) )  ## makes sure x are unit vectors
    apostasi <- tcrossprod( x )
    apostasi <- as.matrix(apostasi)
    diag(apostasi) <- 1
    apostasi[ apostasi >= 1 ] <- 1
    apostasi <- acos(apostasi)
    diag(apostasi) <- 0

    ## The k-NN algorithm is calculated R times. For every repetition a
    ## test sample is chosen and its observations are classified
    for (vim in 1:M) {

      ytest <- as.matrix( y[mat[, vim], ] )  ## test set dependent vars
      ytrain <- as.matrix( y[-mat[, vim], ] )  ## train set dependent vars

      aba <- as.vector( mat[, vim] )
      aba <- aba[aba > 0]
      apo <- apostasi[aba, -aba]
      est <- matrix(nrow = rmat, ncol = d)

      for ( l in 1:c(A - 1) ) {
        k <- l + 1

        if (estim == "arithmetic") {
          for (i in 1:rmat) {
            xa <- cbind(ina, apo[i, ])
            qan <- xa[order(xa[, 2]), ]
            a <- qan[1:k, 1]
            yb <- as.matrix( y[a, ] )
            est[i, ] <- colMeans( yb )
          }

        } else if (estim == "harmonic") {
          for (i in 1:rmat) {
            xa <- cbind(ina, apo[i, ])
            qan <- xa[order(xa[, 2]), ]
            a <- qan[1:k, 1]
            yb <- as.matrix( y[a, ] )
            est[i, ] <- k / colSums( yb )
          }
        }

        if (res == "spher") {
          est <- est / sqrt( rowSums(est^2) )
          ytest <- ytest / sqrt( rowSums(ytest^2) )
          per[vim, l] <- 1 - sum( diag( crossprod(est, ytest) ) ) / rmat
        } else  {
          per[vim, l] <- sum( (est - ytest)^2 ) / rmat
        }
      }
    }

    mspe <- colMeans(per)
    bias <- per[ ,which.min(mspe)] - apply(per, 1, min)  ## TT estimate of bias
    estb <- mean( bias )  ## TT estimate of bias
    performance <- c( 1 - min(mspe) + estb, estb)
    mspe <- 1 - colMeans(per)
    runtime <- proc.time() - runtime

  } else {

    if (ncores == 1) {

      runtime <- proc.time()

      for (vim in 1:M) {
        ytest <- as.matrix( y[mat[, vim], ] )  ## test set dependent vars
        ytrain <- as.matrix( y[-mat[, vim], ] )  ## train set dependent vars
        xtrain <- as.matrix( x[-mat[, vim], ] )  ## train set independent vars
        xtest <- as.matrix( x[mat[, vim], ] )  ## test set independent vars

        for ( l in 1:c(A - 1) ) {
          knn <- l + 1
          est <- knn.reg(xtest, ytrain, xtrain, knn, res = res, type = type, estim = estim)
          per[vim, l] <- sum( (ytest - est)^2 ) / rmat
        }
      }
      runtime <- proc.time() - runtime

    } else {

      runtime <- proc.time()
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      pe <- numeric(A - 1)

      per <- foreach(i = 1:M, .combine = rbind, .export = "knn.reg") %dopar% {
        ytest <- as.matrix( y[mat[, i], ] )  ## test set dependent vars
        ytrain <- as.matrix( y[-mat[, i], ] )  ## train set dependent vars
        xtrain <- as.matrix( x[-mat[, i], ] )  ## train set independent vars
        xtest <- as.matrix( x[mat[, i], ] )  ## test set independent vars
        for ( l in 1:c(A - 1) ) {
          knn <- l + 1
          est <- knn.reg(xtest, ytrain, xtrain, knn, res = res, type = type, estim = estim)
          pe[l] <- sum( (ytest - est)^2 ) / rmat
        }
        return(pe)
      }
      stopCluster(cl)
      runtime <- proc.time() - runtime

    }

    mspe <- colMeans(per)
    bias <- per[ ,which.min(mspe)] - apply(per, 1, min)  ## TT estimate of bias
    estb <- mean( bias )  ## TT estimate of bias
    names(mspe) <- paste("k=", 2:A, sep = "")
    performance <- c( min(mspe) + estb, estb)
    names(performance) <- c("mspe", "estimated bias")

  }

  if (graph == TRUE) {
    plot(2:c(length(mspe) + 1), mspe, xlab = "Nearest neighbours",
    ylab = "MSPE", type = "b")
  }

  list(crit = mspe, best_k = which.min(mspe) + 1, performance = performance, runtime = runtime)
}





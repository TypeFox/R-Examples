################################
#### Multivariate or univariate regression with compositional data
#### in the covariates side using the alpha-transformation
#### Tuning the values of a and the number of principal components
#### Tsagris Michail 2/2016
#### mtsagris@yahoo.gr
#### References:Tsagris M. (2015)
#### Regression analysis with compositional data containing zero values
#### Chilean journal of statistics 6(2): 47-57
################################

alfaridge.tune <- function(y, x, M = 10, a = seq(-1, 1, by = 0.1), lambda = seq(0, 2, by = 0.1),
                           mat = NULL, ncores = 1, graph = TRUE, col.nu = 15) {

  x <- as.matrix(x)
  x <- x/ rowSums(x)
  d <- ncol(x) - 1
  if ( min(x) == 0 )  a <- a[a>0]  ## checks for zero values in the data.
  da <- length(a)
  n <- nrow(x)

  if ( is.null(mat) ) {
    nu <- sample(1:n, min( n, round(n / M) * M ) )
    ## It may be the case this new nu is not exactly the same
    ## as the one specified by the user
    ## to a matrix a warning message should appear
    options(warn = -1)
    mat <- matrix( nu, ncol = M ) # if the length of nu does not fit
  } else  mat <- mat

  M <- ncol(mat)

  mspe2 <- array( dim = c( M, length(lambda), da ) )

  if (ncores <= 1 ) {

    tac <- proc.time()

    for ( i in 1:da ) {
     z <- alfa(x, a[i])$aff
     mod <- ridge.tune(y, z, M = M, lambda = lambda, mat = mat, ncores = 1, graph = FALSE)
     mspe2[, , i] <- mod$msp
    }
    runtime <- proc.time() - tac

  } else {

    tac <- proc.time()
    ## dimensions of the matrix val a warning message should appear
    ## but with options(warn = -1) you will not see it
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    ms <- numeric( M * length(lambda) )
    ww <- foreach(i = 1:da, .combine = cbind, .export = c("ridge.tune",
       "alfa", "helm") ) %dopar% {
      z <- alfa(x, a[i])$aff
      mod <- ridge.tune(y, z, M = M, lambda = lambda, mat = mat, ncores = 1, graph = FALSE)
      ms[i] <- as.vector(mod$msp)
    }
    for (i in 1:da) {
      mspe2[, , i] <- matrix(ww[, i], nrow = M)
    }
    runtime <- proc.time() - tac
  }


  dimnames(mspe2) <- list(folds = 1:M, lambda = lambda, a = a)
  mspe <- array( dim = c(da, length(lambda), M) )
  for (i in 1:M)  mspe[, , i] <- t( mspe2[i, , 1:da] )
  dimnames(mspe) <- list(a = a, lambda = lambda, folds = 1:M )
  mean.mspe <- apply(mspe, 1:2, mean)
  best.par <- ( which(mean.mspe == min(mean.mspe), arr.ind = TRUE)[1, ] )
  opt.mspe <- mean.mspe[ best.par[1], best.par[2] ]
  estb <- mspe[ best.par[1], best.par[2], 1:M ] - apply(mspe, 3, min)
  bias <- mean(estb)
  rownames(mean.mspe) = a   ;  colnames(mspe) = lambda

  if (graph == TRUE) {
    filled.contour( a, lambda, mean.mspe,
                   xlab = expression( paste(alpha, " values") ),
                   ylab = expression( paste(lambda, " values") ) )
  }

  best.par <- c( a[ best.par[1] ], best.par[2] )
  names(best.par) <- c("alpha", "PC")
  performance <- c(opt.mspe, bias)
  names(performance) <- c("bias corrected mspe", "estimated bias")

  runtime <- proc.time()- tac

  list(mspe = mean.mspe, best.par = best.par, performance = performance, runtime = runtime)
}

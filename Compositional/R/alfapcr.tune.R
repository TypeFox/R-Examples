################################
#### Multivariate or univariate regression with compositional data
#### in the covariates side using the alpha-transformation
#### Tuning the values of a and the number of principal components
#### Tsagris Michail 8/2015
#### mtsagris@yahoo.gr
#### References:Tsagris M. (2015)
#### Regression analysis with compositional data containing zero values
#### Chilean journal of statistics 6(2): 47-57
################################

alfapcr.tune <- function(y, x, M = 10, maxk = 50, a = seq(-1, 1, by = 0.1),
  mat = NULL, ncores = 1, graph = TRUE, col.nu = 15) {
  ## oiko can be either "normal", "binomial" or "poisson"
  ## depending on the type of the independent variable
  ## "normal" is set by default

  x <- as.matrix(x)
  x <- x/ rowSums(x)
  n <- nrow(x)
  d <- ncol(x) - 1
  if ( min(x) == 0 )  a <- a[a>0]  ## checks for zero values in the data.
  da <- length(a)

  if ( is.null(mat) ) {
    nu <- sample(1:n, min( n, round(n / M) * M ) )
    ## It may be the case this new nu is not exactly the same
    ## as the one specified by the user
    ## to a matrix a warning message should appear
    options(warn = -1)
    mat <- matrix( nu, ncol = M )
  } else  mat <- mat

  M <- ncol(mat)
  mspe2 <- array( dim = c( M, d, da) )

  if ( length(unique(y)) == 2 ) {
    oiko <- "binomial"
  } else if ( sum(y) - round(y) == 0 ) {
    oiko <- "poisson"
  } else oiko <- "normal"

  if (oiko == 'normal') {

    if ( ncores <=1 ) {
      tic <- proc.time()
      for ( i in 1:da ) {
        z <- alfa(x, a[i])$aff
        mod <- pcr.tune(y, z, M = M, maxk = maxk, mat = mat,
        ncores = 1, graph = FALSE)
        mspe2[, , i] <- mod$msp
      }
      toc <- proc.time() - tic

    } else {

      tac <- proc.time()
      ## dimensions of the matrix val a warning message should appear
      ## but with options(warn = -1) you will not see it
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      ms <- numeric( M * maxk)
      ww <- foreach(i = 1:da, .combine = cbind, .export = c("pcr.tune",
         "alfa", "helm") ) %dopar% {
        z <- alfa(x, a[i])$aff
        mod <- pcr.tune(y, z, M = M, maxk = maxk, mat = mat, ncores = 1, graph = FALSE)
        ms[i] <- as.vector(mod$msp)
      }
      for (i in 1:da) {
        mspe2[, , i] <- matrix(ww[, i], nrow = M)
      }
      runtime <- proc.time() - tac

    }

  } else {

    if (ncores <= 1) {
      tac <- proc.time()

      for ( i in 1:da ) {
        z <- alfa(x, a[i])$aff
        mod <- glmpcr.tune(y, z, M = M, maxk = maxk,
        mat = mat, ncores = ncores, graph = FALSE)
        mspe2[, , i] <- mod$msp
      }
      runtime <- proc.time() - tac

    } else {
      tac <- proc.time()
      ## dimensions of the matrix val a warning message should appear
      ## but with options(warn = -1) you will not see it
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      ms <- numeric( M * maxk)
      ww <- foreach(i = 1:da, .combine = cbind, .export = c("glmpcr.tune",
         "alfa", "helm") ) %dopar% {
        z <- alfa(x, a[i])$aff
        mod <- glmpcr.tune(y, z, M = M, maxk = maxk, mat = mat, ncores = 1, graph = FALSE)
        ms[i] <- as.vector(mod$msp)
                                                            }
      for (i in 1:da) {
        mspe2[, , i] <- matrix(ww[, i], nrow = M)
      }
      runtime <- proc.time() - tac

    }

  }

  dimnames(mspe2) <- list(folds = 1:M, PC = paste("PC", 1:d, sep = ""), a = a)
  mspe <- array( dim = c(da, d, M) )
  for (i in 1:M)  mspe[, , i] <- t( mspe2[i, , 1:da] )
  dimnames(mspe) <- list(a = a, PC = paste("PC", 1:d, sep = ""), folds = 1:M )
  mean.mspe <- apply(mspe, 1:2, mean)
  best.par <- ( which(mean.mspe == min(mean.mspe), arr.ind = TRUE)[1, ] )
  opt.mspe <- mean.mspe[ best.par[1], best.par[2] ]
  estb <- mspe[ best.par[1], best.par[2], 1:M ] - apply(mspe, 3, min)
  bias <- mean(estb)
  rownames(mean.mspe) = a   ;  colnames(mspe) = paste("PC", 1:d, sep = "")

  if (graph == TRUE) {
    filled.contour(a, 1:d, mean.mspe, col = terrain.colors(col.nu),
     xlab = expression( paste(alpha, " values") ), ylab = "Number of PCs")
  }

  best.par <- c( a[ best.par[1] ], best.par[2] )
  names(best.par) <- c("alpha", "PC")
  performance <- c(opt.mspe, bias)
  names(performance) <- c("bias corrected mspe", "estimated bias")
  list(mspe = mean.mspe, best.par = best.par, performance = performance, runtime = toc)
}


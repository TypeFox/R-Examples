################################
#### Tuning the alfa in alfa-regression via K-fold cross validation
#### The bias corrected performance is returned using the
#### Tibshirani and Tibshirani method
#### Tsagris Michail 11/2015
#### mtsagris@yahoo.gr
#### References: Tsagris Michail (2015)
#### Regression analysis with compositional data containing zero values
#### Chilean Journal of Statistics, 6(2): 47-57
#### Tibshirani and Tibshirani (2009),
#### A bias correction for the minimum error rate in cross-validation
#### The Annals of Applied Statistics, 3(1):822-829
################################

alfareg.tune <- function(y, x, a = seq(0.1, 1, by = 0.1), K = 10, mat = NULL,
                         nc = 1, graph = FALSE) {
  ## y is the compositional data (dependent variable)
  ## x is the independent variables
  ## a is a range of values of alpha
  ## K is the number of folds for the K-fold cross validation
  ## nc is how many cores you want to use, default value is 2
  if ( min(y) == 0 )  a <- a[a>0]
  la <- length(a)
  n <- nrow(y)
  x <- as.matrix(x)
  y <- as.matrix(y)
  y <- y /rowSums(y)

  if ( is.null(mat) ) {
    nu <- sample(1:n, min( n, round(n / K) * K ) )
    ## It may be the case this new nu is not exactly the same
    ## as the one specified by the user
    ## to a matrix a warning message should appear
    options(warn = -1)
    mat <- matrix( nu, ncol = K ) # if the length of nu does not fit
  } else  mat <- mat

  K <- ncol(mat)
  rmat <- nrow(mat)

  if (nc == 1) {
    apa <- proc.time()
    kula <- matrix(nrow = K, ncol = la)
    for (j in 1:la) {
      ytr <- alfa(y, a[j])$aff
      for (i in 1:K) {
        xu <- x[ mat[, i], ]
        yu <- y[ mat[, i], ]
        xa <- x[ -mat[, i], ]
        yb <- ytr[ -mat[, i], ]
        mod <- alfa.reg(yu, xa, a[j], xnew = xu, yb = yb)
        yest <- mod$est
        kula[i, j] <- 2 * sum(yu * log(yu / yest), na.rm = TRUE)
      }
    }

    kl <- colMeans(kula)
    opt <- a[ which.min(kl) ]
    val <- which.min(kl)
    per <- min(kl)
    pera <- apply(kula, 1, min)
    bias <- mean( kula[, val] - pera )
    apa <- proc.time() - apa

  } else {
    apa <- proc.time()
    options(warn = -1)
    val <- matrix(a, ncol = nc) ## if the length of a is not equal to the
    ## dimensions of the matrix val a warning message should appear
    ## but with options(warn = -1) you will not see it
    cl <- makePSOCKcluster(nc)
    registerDoParallel(cl)
    kula <- foreach(j = 1:nc, .combine = cbind,
                    .export = c("alfa.reg", "alfa", "helm", "comp.reg", "multivreg") ) %dopar% {
                      ba <- val[, j]
                      ww <- matrix(nrow = K, ncol = length(ba) )
                      for ( l in 1:length(ba) ) {
                        ytr <- alfa(y, ba[l])$aff
                        for (i in 1:K) {
                          xu <- x[ mat[, i], ]
                          yu <- y[ mat[, i], ]
                          xa <- x[ -mat[, i], ]
                          yb <- ytr[ -mat[, i], ]
                          mod <- alfa.reg(yu, xa, ba[l], xnew = xu, yb = yb)
                          yest <- mod$est
                          ww[i, l] <- 2 * sum(yu * log(yu / yest), na.rm = T)
                        }
                      }
                      return(ww)
                    }

    stopCluster(cl)
    kula <- kula[, 1:la]
    kl <- colMeans(kula)
    opt <- a[ which.min(kl) ]
    val <- which.min(kl)
    per <- min(kl)
    pera <- apply(kula, 1, min)
    bias <- mean( kula[, val] - pera )
    apa <- proc.time() - apa
  }

  if (graph == TRUE) {
    plot(a, kula[1, ], type = 'l', ylim = c( min(kula), max(kula) ), xlab = expression(alpha),
         ylab = 'Twice the Kullback Leibler divergence')
    for (i in 2:K)  lines(a, kula[i, ])
    lines(a, kl, col = 2, lty = 2, lwd = 2)
  }

  list(runtime = apa, kl = kl, opt = opt, value = per + bias, bias = bias)
}

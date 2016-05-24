rda.tune <- function(x, ina, M = 10, gam = seq(0, 1, by = 0.1),
                     del = seq(0, 1, by = 0.1), ncores = 1, mat = NULL) {

  ## x contains the data
  ## gam is between pooled covariance and diagonal
  ## gam*Spooled+(1-gam)*diagonal
  ## del is between QDA and LDA
  ## del*QDa+(1-del)*LDA
  ## if ncores==1, then 1 processor is used, otherwise more are
  ## used (parallel computing)
  ## if a matrix with folds is supplied in mat the results will
  ## always be the same. Leave it NULL otherwise

  x <- as.matrix(x)
  ina <- as.numeric(ina)
  n <- nrow(x)  ## total sample size
  num <- 1:n
  nc <- max(ina) ## number of groups
  D <- ncol(x)  ## number of variables
  Ska <- array( dim = c(D, D, nc) )
  ng <- as.vector( table(ina) )
  ci <- log(ng / n)
  con <-  - D / 2 * log(2 * pi)  ## constant part
  sk <- array( dim = c(D, D, nc) )
  lg <- length(gam)  ;  ld <- length(del)

  if ( is.null(mat) ) {
    nu <- sample(1:n, min( n, round(n / M) * M ) )
    ## It may be the case this new nu is not exactly the same
    ## as the one specified by the user
    ## to a matrix a warning message should appear
    options(warn = -1)
    mat <- matrix( nu, ncol = M ) # if the length of nu does not fit
  } else  mat <- mat

  ## mat contains the positions of the test set
  ## this is stored but not showed in the end
  ## the user can access it though by running
  ## the commands outside this function
  rmat <- nrow(mat)
  M <- ncol(mat)
  gr <- matrix(nrow = rmat, ncol = nc)
  msp <- array(dim = c(lg, ld, M) )

  if (ncores > 1) {
    runtime <- proc.time()
    group <- matrix(nrow = length(gam), ncol = length(del) )
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)

    ww <- foreach(vim = 1:M, .combine = cbind) %dopar% {

      test <- as.matrix( x[ mat[, vim], ] )  ## test sample
      id <- as.vector( ina[ mat[, vim] ] )  ## groups of test sample
      train <- as.matrix( x[ -mat[, vim], ] )   ## training sample
      ida <- as.vector( ina[ -mat[, vim] ] )   ## groups of training sample

      na <- as.vector( table(ida) )
      na <- rep(na - 1, each = D^2)
      mesi <- aggregate( train, by = list(ida), mean )
      mesi <- as.matrix( mesi[, -1] )

      ## the covariance matrix of each group is now calculated
      for (m in 1:nc)  sk[ , , m] <- cov( train[ida == m, ] )
      s <- na * sk
      Sp <- apply(s, 1:2, sum) / (n - nc)  ## pooled covariance matrix
      sp <- diag( mean( diag( Sp ) ), D )

      for (k1 in 1:length(gam)) {
        for (k2 in 1:length(del)) {
          Sa <- gam[k1] * Sp + (1 - gam[k1]) * sp  ## regularised covariance matrix
          for (j in 1:nc) {
            Ska[, , j] <- del[k2] * sk[, , j] + (1 - del[k2]) * Sa
            gr[, j] <- ci[j] - 0.5 * log( det( Ska[, , j] ) ) -
              0.5 * mahalanobis( test, mesi[j, ], Ska[, , j] )
          }
          gr <- gr + con
          g <- apply(gr, 1, which.max)
          group[k1, k2] <- mean(g == id)
        }
      }
      a <- as.vector( group )
      return(a)
    }
    stopCluster(cl)

    per <- array( dim = c( lg, ld, M ) )
    index <- matrix( 1:c(lg * ld * M), ncol = M )

    for ( i in 1:M ) {
      per[, , i] <- matrix( ww[, i], nrow = lg )
    }

    runtime <- proc.time() - runtime

  } else {
    runtime <- proc.time()
    per <- array( dim = c( lg, ld, M ) )

    for (vim in 1:M) {

      test <- as.matrix( x[ mat[, vim], ] )  ## test sample
      id <- as.vector( ina[ mat[, vim] ] )  ## groups of test sample
      train <- as.matrix( x[ -mat[, vim], ] )   ## training sample
      ida <- as.vector( ina[ -mat[, vim] ] )   ## groups of training sample

      na <- as.vector( table(ida) )
      na <- rep(na - 1, each = D^2)
      mesi <- aggregate( train, by = list(ida), mean )
      mesi <- as.matrix( mesi[, -1] )

      ## the covariance matrix of each group is now calculated
      for (m in 1:nc)  sk[ , , m] <- cov( train[ida == m, ] )
      s <- na * sk
      Sp <- apply(s, 1:2, sum) / (n - nc)  ## pooled covariance matrix
      sp <- diag( mean( diag( Sp ) ), D )

      for (k1 in 1:length(gam)) {
        for (k2 in 1:length(del)) {
          Sa <- gam[k1] * Sp + (1 - gam[k1]) * sp  ## regularised covariance matrix
          for (j in 1:nc) {
            Ska[, , j] <- del[k2] * sk[, , j] + (1 - del[k2]) * Sa
            gr[, j] <- ci[j] - 0.5 * log( det( Ska[, , j] ) ) -
              0.5 * mahalanobis( test, mesi[j, ], Ska[, , j] )
          }
          gr <- gr + con
          g <- apply(gr, 1, which.max)
          per[k1, k2, vim] <- mean(g == id)
        }
      }
    }
    runtime <- proc.time() - runtime

  }

  percent <- apply(per, 1:2, mean)
  su <- apply(per, 1:2, sd)
  dimnames(percent) <- dimnames(su) <- list(gamma = gam, delta = del)

  confa <- as.vector( which(percent == max( percent ), arr.ind = TRUE )[1, ] )
  bias <- numeric(M)
  for (i in 1:M) {
    confi <- as.vector( which(per[, , i] == max( per[, , i] ), arr.ind = TRUE )[1, ] )
    bias[i] <- per[ confi[1], confi[2], i] - per[ confa[1], confa[2], i]
  }

  result <- cbind( max(percent) - mean(bias), gam[ confa[1] ], del[ confa[2] ] )
  colnames(result) <- c('optimal', 'best gamma', 'best delta')

  list(per = per, percent = percent, se = su, result = result, runtime = runtime)
}

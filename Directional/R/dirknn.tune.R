################################
#### Discriminant analysis for directional data
#### using the k-NN alorithm, tuning the k neighbours
#### Tsagris Michail 01/2016
#### mtsagris@yahoo.gr
################################

dirknn.tune <- function(z, M = 10, A = 5, ina, type = "S",
                        mesos = TRUE, mat = NULL) {
  ## x is the matrix containing the data
  ## M is the number of folds, set to 10 by default
  ## A is the maximum number of neighbours to use
  ## ina indicates the groups, numerical variable
  ## type is either 'S' or 'NS'. Should the standard k-NN be use or not
  ## if mesos is TRUE, then the arithmetic mean distange of the k nearest
  ## points will be used.
  ## If not, then the harmonic mean will be used. Both of these apply for
  ## the non-standard algorithm, that is when type='NS'

  runtime <- proc.time()
  z <- as.matrix(z)  ## makes sure the x is a matrix
  z <- z / sqrt( rowSums(z^2) )  ## makes sure the the data are unit vectors
  n <- nrow(z)  ## sample size
  ina <- as.numeric(ina)
  if ( A >= min(table(ina)) )  A <- min(table(ina)) - 3  ## The maximum
  ## number  of nearest neighbours to use
  ina <- as.numeric(ina) ## makes sure ina is numeric
  ng <- max(ina)  ## The number of groups

  if ( is.null(mat) ) {
    nu <- sample(1:n, min( n, round(n / M) * M ) )
    ## It may be the case this new nu is not exactly the same
    ## as the one specified by the user
    ## to a matrix a warning message should appear
    options(warn = -1)
    mat <- matrix( nu, ncol = M )
  } else  mat <- mat

  M <- ncol(mat)
  per <- matrix(nrow = M, ncol = A)
  rmat <- nrow(mat)

  dis <- tcrossprod( z )
  diag(dis) <- 1
  dis[ dis > 1 ] <- 1
  dis <- acos(dis)

  ## The k-NN algorith is calculated M times. For every repetition a
  ## fold is chosen and its observations are classified
  for (vim in 1:M) {

    id <- as.vector( ina[ mat[, vim] ] )  ## groups of test sample
    ina2 <- as.vector( ina[ -mat[, vim] ] )   ## groups of training sample
    aba <- as.vector( mat[, vim] )
    aba <- aba[aba > 0]
    apo <- dis[aba, -aba]
    ta <- matrix(nrow = rmat, ncol = ng)

    if (type == "NS") {
      ## Non Standard algorithm
      for ( j in 1:c(A - 1) ) {
        knn <- j + 1
        for (l in 1:ng) {
          dista <- apo[, ina2 == l]
          dista <- t( apply(dista, 1, sort) )
          if (mesos == TRUE) {
            ta[, l] <- rowMeans( dista[, 1:knn] )
          } else {
            ta[, l] <- knn / rowSums( 1 / dista[, 1:knn] )
          }
        }
        g <- apply(ta, 1, which.min)
        per[vim, j] <- mean(g == id)
      }

    } else {
      ## Standard algorithm
      g <- numeric( A - 1)
      for ( j in 1:c(A - 1) ) {
        knn <- j + 1
        for (k in 1:rmat) {
          xa <- cbind(ina2, apo[k, ])
          qan <- xa[order(xa[, 2]), ]
          sa <- qan[1:knn, 1]
          tab <- table(sa)
          g[k] <- as.integer( names(tab)[ which.max(tab) ] )
        }
        per[vim, j] <- mean(g == id)
      }
    }
  }

  ela <- colMeans(per)
  bias <- per[ , which.max(ela)] - apply(per, 1, max)  ## TT estimate of bias
  estb <- mean( bias )  ## TT estimate of bias
  runtime <- proc.time() - runtime
  names(ela) <- paste("k=", 2:c(A + 1), sep = "")
  plot(2:c(A + 1), ela, type = "b", xlab = "k nearest neighbours",
       pch = 9, ylab = "Estimated percentage of correct classification")
  percent <- c( max(ela) + estb)
  names(percent) <- c("Bias corrected estimated percentage")

  list( per = ela, percent = percent, runtime = runtime )
}

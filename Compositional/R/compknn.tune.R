################################
#### Classification for compositional data using a power transformation
#### Tuning the k-NN algorithm via M-fold cross-validation
#### Tsagris Michail 7/2015
#### References: Tsagris, M. T. (2014).
#### The k-NN algorithm for compositional data: a revised approach with and without zero values present
#### Journal of Data Science, 12(3):519-534
#### mtsagris@yahoo.gr
################################

compknn.tune <- function(x, ina, M = 10, A = 5, type= "S", mesos = TRUE,
                         a = seq(-1, 1, by = 0.1), apostasi = "ESOV", mat = NULL, graph = FALSE) {

  ## x is the matrix containing the data
  ## M is the number of folds, set to 10 by default
  ## A is the maximum number of neighbours to use
  ## actually it goes until A + 1
  ## ina indicates the groups, numerical variable
  ## a is a vector containing the values of the power parameter
  ## type is either 'S' or 'NS'. Should the standard k-NN be use or not
  ## if mesos is TRUE, then the arithmetic mean distange of the k nearest
  ## points will be used.
  ## If not, then the harmonic mean will be used. Both of these apply for
  ## the non-standard algorithm, that is when type='NS'
  ## apostasi is the type of metric used: 'ESOV' or 'taxicab',
  ## 'Ait', 'Hellinger', or 'angular'

  x <- as.matrix(x)  ## makes sure the x is a matrix
  x <- x/rowSums(x)  ## makes sure the the data sum to 1
  n <- nrow(x)  ## sample size
  ina <- as.numeric(ina)
  if ( A >= min(table(ina)) )  A <- min(table(ina)) - 3  ## The maximum
  ## number  of nearest neighbours to use
  ng <- max(ina)  ## The number of groups

  dis <- matrix(numeric(n^2), nrow = n, ncol = n)
  ## The next two functions split the sample into R different test
  ## and training datasets
  ## The test dataset is chosen via stratified or simple random sampling
  ## will be stored in the array called per

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

  ## The algorithm is repated R times and each time the estimated
  ## percentages are stored in the array per.

  if (apostasi == "ESOV" | apostasi == "taxicab") {

    runtime <- proc.time()
    per <- array( dim = c(M, A, length(a)) )

    for (i in 1:length(a)) {

      z <- x^a[i] / rowSums( x^a[i] )  ## The power transformation is applied

      if (apostasi == "ESOV") {
        for ( m1 in 1:c(n - 1) ) {
          for ( m2 in c(m1 + 1):n ) {
            ma <- z[m1, ] + z[m2, ]
            dis[m1, m2] <- sqrt( sum( z[m1, ] * log( 2 * z[m1, ]/ma ) +
                                        z[m2, ] * log( 2 * z[m2, ]/ma ), na.rm = TRUE ) )
          }
        }
        dis <- dis + t(dis)

      } else if (apostasi == "taxicab") {
        dis <- dist(z, method = "manhattan", diag = TRUE, upper = TRUE)
        dis <- as.matrix(dis)
      }

      ## The k-NN algorithm is calculated R times. For every repetition a
      ## test sample is chosen and its observations are classified
      for (vim in 1:M) {

        id <- as.vector( ina[ mat[, vim] ] )  ## groups of test sample
        ina2 <- as.vector( ina[ -mat[, vim] ] )   ## groups of training sample
        aba <- as.vector( mat[, vim] )
        aba <- aba[aba > 0]
        apo <- dis[aba, -aba]
        ta <- matrix(nrow = rmat, ncol = ng)

        if (type == "NS") {
          ## Non Standard algorithm
          for (j in 1:A) {
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
            per[vim, j, i] <- mean(g == id)
          }

        } else if (type == "S") {
          ## Standard algorithm
          for (j in 1:A) {
            g <- numeric(rmat)
            knn <- j + 1
            for (k in 1:rmat) {
              xa <- cbind(ina2, apo[k, ])
              qan <- xa[order(xa[, 2]), ]
              sa <- qan[1:knn, 1]
              tab <- table(sa)
              g[k] <- as.integer(names(tab)[which.max(tab)])
            }
            per[vim, j, i] <- mean(g == id)
          }
        }
      }

    }

    ela <- matrix(nrow = length(a), ncol = A)
    for ( i in 1:length(a) )  ela[i, ] <- colMeans(per[, , i])
    ## The ela matrix contains the averages of the R
    ## repetitions over alpha and k
    colnames(ela) <- paste("k=", 2:c(A + 1), sep = "")
    rownames(ela) <- paste("alpha=", a, sep = "")

    ## The code for the heat plot of the estimated percentages
    if (graph == TRUE) {
      fields::image.plot(a, 2:c(A + 1), ela, col = grey(1:11/11),
                         ylab = "k nearest-neighbours",
                         xlab = expression(paste(alpha, " values")) )
    }

    opt <- max(ela)
    confa <- as.vector( which(ela == opt, arr.ind = TRUE)[1, ] )
    bias <- numeric(M)
    for (i in 1:M) {
      bias[i] <- opt - per[ i, confa[2], confa[1] ]
    }
    bias <- mean(bias)
    performance <- c(opt - bias, bias)
    names(performance) <- c( "rate", "bias" )
    runtime <- proc.time() - runtime
    results <- list( ela = ela, performance = performance,
                     best_a = a[ confa[1] ], best_k = confa[2] + 1, runtime = runtime )


  } else if (apostasi == "Ait" | apostasi == "Hellinger" | apostasi == "angular") {

    runtime <- proc.time()
    per <- matrix(nrow = M, ncol = A)

    if (apostasi == "Ait") {
      xa <- log(x)
      z <- xa - rowMeans( xa )
      dis <- fields::rdist(z)

    } else if (apostasi == "Hellinger") {
      z <- sqrt(x)
      dis <- fields::rdist(z)
      dis <- dis / sqrt(2)

    } else if (apostasi == "angular") {
      z <- sqrt(x)
      dis <- tcrossprod( z )
      diag(dis) <- 1
      dis[ dis > 1 ] <- 1
      dis <- acos(dis)
    }
    diag(dis) <- 0

    for (vim in 1:M) {

      id <- as.vector( ina[ mat[, vim] ] )  ## groups of test sample
      ina2 <- as.vector( ina[ -mat[, vim] ] )   ## groups of training sample
      aba <- as.vector( mat[, vim] )
      aba <- aba[aba > 0]
      apo <- dis[aba, -aba]
      ta <- matrix(nrow = rmat, ncol = ng)

      if (type == "NS") {
        ## Non Standard algorithm
        for (j in 1:A) {
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
          per[vim, j] <- sum(g == id)/nu
        }

      } else {   ## if (type == "S")
        ## Standard algorithm
        for (j in 1:A) {
          knn <- j + 1
          g <- numeric(rmat)
          for (k in 1:rmat) {
            xa <- cbind(ina2, apo[k, ])
            qan <- xa[order(xa[, 2]), ]
            sa <- qan[1:knn, 1]
            tab <- table(sa)
            g[k] <- as.integer(names(tab)[which.max(tab)])
          }
          per[vim, j] <- mean(g == id)
        }
      }
    }

    ela <- colMeans(per)
    opt <- max(ela)
    names(ela) <- paste("k=", 2:c(A + 1), sep = "")
    best_k = which.max(ela) + 1
    bias <- apply(per, 1, max) - per[, best_k]
    bias <- mean(bias)
    performance <- c(opt - bias, bias)
    names(performance) <- c( "rate", "bias" )

    if (graph == TRUE) {
      plot(2:c(A + 1), ela, type = "b", xlab = "k nearest neighbours", pch = 9,
           col = 2, ylab = "Estimated percentage of correct classification")
    }

    runtime <- proc.time() - runtime

    results <- list(ela = ela, performance = performance, best_k = which.max(ela) + 1,
                    runtime = runtime)
  }

  results
}

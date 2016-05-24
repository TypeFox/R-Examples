################################
#### Classification for compositional data using a power transformation
#### The k-NN algorithm
#### Tsagris Michail 7/2015
#### References: Tsagris, M. T. (2014).
#### The k-NN algorithm for compositional data: a revised approach with and without zero values present
#### Journal of Data Science, 12(3):519-534
#### mtsagris@yahoo.gr
################################

comp.knn <- function(xnew, x, ina, a = 1, k = 5, type = "S",
                     apostasi = "ESOV", mesos = TRUE) {
  ## xnew is the new dataset. It can be a single vector or a matrix
  ## x is the matrix containing the data
  ## ina indicates the groups
  ## a is the value of the power parameter
  ## k in the number of nearest neighbours
  ## apostasi is the type of metric used, "ESOV", "taxicab",
  ## "Ait", "Hellinger", or "angualr"
  ## type is either S or NS. Should the standard k-NN be use or not
  ## if mesos is TRUE, then the arithmetic mean distance of the
  ## k nearest points will be used
  ## If not, then the harmonic mean will be used.
  ## Both of these apply for the non-standard,
  ## algorithm, that is when type=NS

  x <- as.matrix(x)  ## makes sure x is a matrix
  x <- x/rowSums(x)  ## makes sure the data sum to 1
  n <- nrow(x)
  ina <- as.numeric(ina)
  xnew <- as.matrix(xnew)
  xnew <- matrix(xnew, ncol = ncol(x)) ## makes sure xnew is a matrix
  xnew <- xnew/rowSums(xnew)  ## make the data sum to 1
  nc <- max(ina)  ## The number of groups
  nu <- nrow(xnew)
  w <- rbind(x, xnew)

  if (apostasi == "ESOV") {
    nz <- nrow(w)
    dis <- matrix( numeric(nz^2), nrow = nz, ncol = nz )
    z <- w^a / rowSums( w^a )  ## The power transformation is applied
    for (m1 in 1:c(nz - 1)) {
      for (m2 in c(m1 + 1):nz) {
        ma <- z[m1, ] + z[m2, ]
        dis[m1, m2] <- sqrt( sum( z[m1, ] * log( 2 * z[m1, ]/ma ) +
                                    z[m2, ] * log( 2 * z[m2, ]/ma ), na.rm = TRUE ) )
      }
    }
    dis <- dis + t(dis)

  } else  if (apostasi == "taxicab") {
    z <- w^a / rowSums( w^a )  ## The power transformation is applied
    dis <- dist(z, method = "manhattan", diag = TRUE, upper = TRUE)
    dis <- as.matrix(dis)

  } else if (apostasi == "Ait") {
    ## this requires non zero data ## be careful
    xa <- log(x)
    z <- xa - rowMeans( xa )
    dis <- fields::rdist(z)

  } else if (apostasi == "Hellinger") {
    z <- sqrt(x)
    dis <- fields::rdist(z)
    dis <- dis /sqrt(2)

  } else if (apostasi == "angular") {
    z <- sqrt(x)
    dis <- tcrossprod( z )
    diag(dis) <- 1
    dis[dis > 1] <- 1
    dis <- acos(dis)
  }

  ta <- matrix(nrow = nu, ncol = nc)
  apo <- matrix(dis[-c(1:n), 1:n], nrow = nu)

  if (type == "NS") {
    ## Non Standard algorithm
    for (m in 1:nc) {
      dista <- apo[, ina == m]
      dista <- t( apply(dista, 1, sort) )
      if (mesos == TRUE) {
        ta[, m] <- rowMeans( dista[, 1:k] )

      } else {
        ta[, m] <- k / rowSums( 1 / dista[, 1:k] )
      }
    }
    g <- apply(ta, 1, which.min)

  } else {   ## if type is "S"
    ## Standard algorithm
    g <- numeric(nu)
    for (l in 1:nu) {
      xa <- cbind(ina, apo[l, ])
      qan <- xa[order(xa[, 2]), ]
      sa <- qan[1:k, 1]
      tab <- table(sa)
      g[l] <- as.integer(names(tab)[which.max(tab)])
    }
  }

  return(g)
}

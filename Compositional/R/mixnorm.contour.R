################################
#### Contour plot of a normal mixture model in S^2
#### Tsagris Michail 5/2015
#### mtsagris@yahoo.gr
################################

mixnorm.contour <- function(x, mod) {
  ## mod is a mixture model containing all the parameters
  ## the type parameter determines whether the additive or the isometric
  ## log-ratio transformation will be used. If type='alr' (the default) the
  ## additive log-ratio transformation is used. If type='ilr', the isometric
  ## log-ratio is used
  x <- as.matrix(x)  ## makes sure x is matrix
  x <- x/rowSums(x)  ## make sure x compositional data
  prob <- mod$prob  ## mixing probabilitiy of each cluster
  mu <- mod$mu
  su <- mod$su
  type <- mod$type  ## the type of the log-ratio transformation, either "alr" or "ilr"
  g <- length(mod$prob)  ## how many clusters are there
  n <- 100  ## n is the number of points of each axis used
  x1 <- seq(0.001, 0.999, length = n)
  x2 <- seq(0.001, sqrt(3)/2 - 1e-04, length = n)
  mat <- matrix(nrow = n, ncol = n)
  ha <- t( helm (3) )
  for (i in 1:c(n/2)) {
    for (j in 1:n) {
      if (x2[j] < sqrt(3) * x1[i]) {
        ## This checks if the point will lie inside the triangle
        ## The next 4 lines calculate the composition
        w3 <- (2 * x2[j]) / sqrt(3)
        w2 <- x1[i] - x2[j] / sqrt(3)
        w1 <- 1 - w2 - w3
        w <- c(w1, w2, w3)
        if (type == "alr") y <- log( w[-3]/w[3] )  ## alr transformation
        if (type == "ilr") {  ## isometric log-ratio transformation
          y <- log(w) - mean(log(w))
          y <- as.vector( y %*% ha )
        }
        ta <- numeric(g)
        for (k in 1:g) {
          ta[k] <- -0.5 * log(det(2 * pi * su[, , k])) -
          0.5 * mahalanobis(y, mu[k, ], su[, , k])
        }
        can <- sum( prob * exp(ta) )
        if (abs(can) < Inf) {
          mat[i, j] <- can
        } else  mat[i, j] <- NA
      }
    }
  }

  for (i in c(n/2 + 1):n) {
    for (j in 1:n) {
      ## This checks if the point will lie inside the triangle
      if (x2[j] < sqrt(3) - sqrt(3) * x1[i]) {
        ## The next 4 lines calculate the composition
        w3 <- (2 * x2[j]) / sqrt(3)
        w2 <- x1[i] - x2[j] / sqrt(3)
        w1 <- 1 - w2 - w3
        w <- c(w1, w2, w3)
        if (type == "alr") y <- log( w[-3]/w[3] )  ## alr transformation
        if (type == "ilr") {  ## isometric log-ratio transformation
          y <- log(w) - mean(log(w))
          y <- as.vector( y %*% ha )
        }
        ta <- numeric(g)
        for (k in 1:g) {
          ta[k] <- -0.5 * log(det(2 * pi * su[, , k])) -
          0.5 * mahalanobis(y, mu[k, ], su[, , k])
        }
        can <- sum( prob * exp(ta) )
        if (abs(can) < Inf) {
          mat[i, j] <- can
        } else  mat[i, j] <- NA
      }
    }
  }

  contour(x1, x2, mat, col = 3)
  b1 <- c(1/2, 0, 1, 1/2)
  b2 <- c(sqrt(3)/2, 0, 0, sqrt(3)/2)
  b <- cbind(b1, b2)
  points(b[ , 1], b[ , 2] , type = "l", xlab = " ", ylab = " ")
  proj <- matrix(c(0, 1, 1/2, 0, 0, sqrt(3)/2), ncol = 2)
  xa <- x %*% proj
  points(xa[, 1], xa[, 2])
}

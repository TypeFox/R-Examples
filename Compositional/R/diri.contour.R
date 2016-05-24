################################
#### Contour plot of the Dirichlet distribution in S^2
#### Tsagris Michail 1/2013
#### mtsagris@yahoo.gr
################################

diri.contour <- function(a, n = 100, x = NULL) {
  ## a are the estimated Dirichlet parameters
  ## n shows the number of points at which the density is calculated
  ## so, n^2 points are used.
  ## x should be a 3-part compositional data or NULL for no data
  x1 <- seq(0.001, 0.999, length = n)  ## coordinates of x
  x2 <- seq(0.001, sqrt(3)/2 - 1e-04, length = n)  ## coordinates of y
  mat <- matrix(nrow = n, ncol = n)
  beta <- prod( gamma(a)) / gamma(sum(a) )  ## beta function
  for (i in 1:c(n/2)) {
    for (j in 1:n) {
      if (x2[j] < sqrt(3) * x1[i]) {
        ## This checks if the point will lie inside the triangle
        ## the next three lines invert the points which lie inside
        ## the triangle back into the composition in S^2
        w3 <- (2 * x2[j])/sqrt(3)
        w2 <- x1[i] - x2[j]/sqrt(3)
        w1 <- 1 - w2 - w3
        w <- c(w1, w2, w3)
        can <- (1 / beta) * prod(w^(a - 1))
        if (abs(can) < Inf)  mat[i, j] <- can  else  mat[i, j] <- NA
      } else  mat[i, j] <- NA
    }
  }
  for (i in c(n/2 + 1):n) {
    for (j in 1:n) {
      ## This checks if the point will lie inside the triangle
      if (x2[j] < sqrt(3) - sqrt(3) * x1[i]) {
        ## the next three lines invert the points which lie inside
        ## the triangle back into the composition in S^2
        w3 <- (2 * x2[j])/sqrt(3)
        w2 <- x1[i] - x2[j]/sqrt(3)
        w1 <- 1 - w2 - w3
        w <- round(c(w1, w2, w3), 6)
        can <- (1 / beta) * prod(w^(a - 1))
        if (abs(can) < Inf)  mat[i, j] <- can  else  mat[i, j] <- NA
      } else  mat[i, j] <- NA
    }
  }
  contour(x1, x2, mat, col = 3)  ## contour plots
  b1 <- c(1/2, 0, 1, 1/2)
  b2 <- c(sqrt(3)/2, 0, 0, sqrt(3)/2)
  b <- cbind(b1, b2)
  ## the next line draws the triangle in the two dimensions
  points(b[, 1], b[, 2], type = "l", xlab = " ", ylab = " ")
  if ( !is.null(x) ) {
    x <- as.matrix(x)
    x <- x / rowSums(x)
    proj <- matrix(c(0, 1, 1/2, 0, 0, sqrt(3)/2), ncol = 2)
    x <- as.matrix(x)
    x <- x/rowSums(x)
    xa <- x %*% proj
    points(xa[, 1], xa[, 2])
  }
}

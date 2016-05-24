################################
#### Contour plots of the bivariate t distribution in S^2
#### Tsagris Michail 2/2013
#### mtsagris@yahoo.gr
################################

bivt.contour <- function(x, type = 'alr', n = 100, appear = TRUE) {
  ## x is the compositional data
  ## the type parameter determines whether the additive
  ## or the isometric log-ratio transformation will be used.
  ## If type='alr' (the default) the additive
  ## log-ratio transformation is used.
  ## If type='ilr', the isometric log-ratio is used
  ## n is the number of points of each axis used
  x <- as.matrix(x)
  x <- x / rowSums(x)
  if (type == 'alr') {
    y <- log( x[, -3] / x[, 3] ) ## additive log-ratio transformation
  } else {
    y <- log(x) - rowMeans( log(x) )
    y <- as.matrix( y %*% t( helm(3) ) )
  }
  mod <- multivt(y)
  m <- mod$center
  s <- mod$scatter
  v <- mod$df
  p <- 2
  x1 <- seq(0.001, 0.999, length = n)
  x2 <- seq(0.001, sqrt(3)/2 - 0.0001, length = n)
  mat <- matrix(nrow = n, ncol = n)
  st <- solve(s)
  for (i in 1:c(n/2) ) {
   for (j in 1:n) {
     if (x2[j] < sqrt(3) * x1[i]) { ## This checks if the point lies
     ## inside the triangle
   ## The next 4 lines calculate the composition
       w3 <- (2 * x2[j]) / sqrt(3)
       w2 <- x1[i] - x2[j] / sqrt(3)
       w1 <- 1 - w2 - w3
       w <- c(w1, w2, w3)
       if (type == 'alr') {
	     y <- log(w[-3] / w[3]) ## additive log-ratio transformation
       } else {
         y <- log(w) - mean(log(w))
         y <- as.vector( y %*% t(helm(3)) )
       }  ## isometric log-ratio transformation
       ca <- lgamma( (v + p)/2 ) - lgamma(v/2) - 0.5 * log( det(pi * v * s) )-
       0.5 * (v + p) * ( log( 1 + ( c(y[1] - m[1], y[2] - m[2]) %*%
       st %*% c( y[1] - m[1], y[2] - m[2]) )/v ) )
       can <- exp(ca)
  	   if (abs(can) < Inf) {
	     mat[i, j] <- can
	   } else mat[i, j] <- NA
	 }
   }
  }
  for (i in c(n/2 + 1):n) {
   for (j in  1:n) {
  ## This checks if the point will lie inside the triangle
    if (x2[j] < sqrt(3) - sqrt(3) * x1[i]) {
  ## The next 4 lines calculate the composition
     w3 <- (2 * x2[j]) / sqrt(3)
     w2 <- x1[i] - x2[j] / sqrt(3)
     w1 <- 1 - w2 - w3
     w <- c(w1, w2, w3)
     if (type == 'alr') {
	   y <- log(w[-3]/w[3]) ## additive log-ratio transformation
     } else  {
       y <- log(w) - mean(log(w))
       y <- as.vector( y %*% t(helm(3)) )
     }  ## isometric log-ratio transformation
     ca <- lgamma((v + p)/2) - lgamma(v/2) - 0.5 * log( det(pi * v * s) ) -
     0.5 * (v + p) * (log( 1 + ( c( y[1] - m[1], y[2] - m[2] ) %*%
     st %*% c( y[1] - m[1], y[2] - m[2]) )/v ) )
     can <- exp(ca)
     if (abs(can) < Inf) {
	  mat[i, j] <- can
	 } else mat[i, j] <- NA
	}
   }
  }
  contour(x1, x2, mat, nlevels = 7, col = 3, pty = "s", xaxt = "n",
  yaxt = "n", bty = "n")
  b1 <- c(1/2, 0, 1, 1/2)
  b2 <- c(sqrt(3)/2, 0, 0, sqrt(3)/2)
  b <- cbind(b1, b2)
  points(b[, 1], b[, 2], type = "l", xlab = " ", ylab = " ")
  if (appear == TRUE){
   nam <- colnames(x)
   text(b[1, 1], b[1, 2] + 0.02, nam[3], cex = 1)
   text(b[2:3, 1], b[2:3, 2] - 0.02, nam[1:2], cex = 1)
   proj <- matrix(c(0, 1, 1/2, 0, 0, sqrt(3)/2), ncol = 2)
   x <- as.matrix(x) ; x <- x/rowSums(x)
   xa <- x %*% proj
   points(xa[, 1], xa[, 2])
  }
}

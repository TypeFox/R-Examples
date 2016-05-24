################################
#### Ternary plot
#### Tsagris Michail 5/2012
#### mtsagris@yahoo.gr
################################

ternary <- function(x, means = TRUE, pca = FALSE) {
  ## x contains the composiitonal data
  ## if means==TRUE it will plot the arithmetic and the
  ## closed geometric mean
  ## if pca==TRUE it will plot the first principal component
  x <- as.matrix(x)  ## makers sure x is a matrix
  x <- x / rowSums(x)  ## makes sure x is compositional data
  if ( !is.null( colnames(x) ) ) {
    nam <- colnames(x)
  } else nam <- paste("X", 1:3, sep = " ")
  n <- nrow(x)
  ina <- numeric(n) + 1
  ## m1 is the closed geometric mean
  g1 <- colMeans( log(x[, -1] / x[, 1]) )
  g2 <- c( 1, exp(g1) )
  m1 <- g2 / sum(g2)
  ## m2 is the simple arithmetic mean
  m2 <- colMeans(x)
  x <- rbind(x, m1, m2)
  ## the next code checks for zeros
  ina[ rowSums(x == 0) == 1 ] <- 3
  b1 <- c(1/2, 0, 1, 1/2)
  b2 <- c(sqrt(3)/2, 0, 0, sqrt(3)/2)
  b <- cbind(b1, b2)
  plot(b[, 1], b[, 2], type = "l", xlab = " ", ylab = " ", pty = "s",
  xaxt = "n", yaxt = "n", bty = "n")
  proj <- matrix(c(0, 1, 1/2, 0, 0, sqrt(3)/2), ncol = 2)
  d <- x %*% proj
  points( d[1:n, 1], d[1:n, 2], col = ina )
  text( b[1, 1], b[1, 2] + 0.02, nam[3], cex = 1 )
  text( b[2:3, 1], b[2:3, 2] - 0.02, nam[1:2], cex = 1 )
  if (means == TRUE) {
    ## should the means appear in the plot?
    points( d[c(n + 1), 1], d[c(n + 1), 2], pch = 2, col = 2 )
    points( d[c(n + 2), 1], d[c(n + 2), 2], pch = 3, col = 3 )
    legend(0.57, 0.9, c("closed geometric mean"," arithmetic mean"),
    pch = c(2, 3), col = c(2, 3), bg = 'gray90')
  }
  if (pca == TRUE  &  min(x) > 0 ) {
    ## should the first principal component appear?
    z <- log( x[1:n, ] ) - rowMeans( log(x[1:n, ]) ) ## clr transformation
    m <- colMeans(z)  ## mean vector in the clr space
    a <- eigen( cov(z) )$vectors[, 1] + m  ## move the unit vector a bit
    sc <- z %*% a
    lam <- seq( min(sc) - 1.5, max(sc) + 1.5, length = n )
    x1 <- cbind( a[1] * lam, a[2] * lam, a[3] * lam) + cbind( m[1] * (1 - lam),
    m[2] * (1 - lam), m[3] * (1 - lam) )
    wa1 <- exp(x1) / rowSums( exp(x1) )  ## first principal component in S^2
    wa <- wa1 %*% proj
    lines(wa, lwd = 2, lty = 2)
  }
  mu <- rbind(m1, m2)
  rownames(mu) <- c("closed geometric", "arithmetic mean")
  colnames(mu) <- nam
  mu
}

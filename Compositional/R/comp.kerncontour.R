################################
#### Contour plot of the kernel density estimate in S^2
#### Tsagris Michail 2/2015
#### mtsagris@yahoo.gr
################################

comp.kerncontour <- function(x, type = "alr", n = 100) {
  ## x contains the compositional data
  ## type determines which log-ratio transformation will be used.
  ## If type='alr' (the default) the additive
  ## log-ratio transformation is used.
  ## If type='ilr', the isometric log-ratio is used
  ## n is the number of points of each axis used

  x <- as.matrix(x)
  x <- x/rowSums(x)  ## makes sure x is a matrix with compositional data
  nu <- nrow(x)  ## sample size
  if (type == "alr")  z <- log(x[, -3]/x[, 3])  ## alr transformation
  if (type == "ilr") {  ## isometric log-ratio transformation
      z <- log(x) - rowMeans(log(x))
      z <- z %*% t(helm(3))
  }

  hopt <- mkde.tune(z)$hopt
  con <- hopt^2
  hopt <- diag( rep(hopt, 2) )
  ts <- solve(hopt^2)
  x1 <- seq(0.001, 0.999, length = n)
  x2 <- seq(0.001, sqrt(3)/2 - 1e-04, length = n)
  mat <- matrix(nrow = n, ncol = n)

  for (i in 1:c(n/2)) {
    for (j in 1:n) {
      if (x2[j] < sqrt(3) * x1[i]) {
        ## This checks if the point will lie inside the triangle
        ## The next 4 lines calculate the composition
        w3 <- (2 * x2[j])/sqrt(3)
        w2 <- x1[i] - x2[j]/sqrt(3)
        w1 <- 1 - w2 - w3
        w <- c(w1, w2, w3)
        if (type == "alr")  y <- log(w[-3]/w[3])  ## alr transformation
        if (type == "ilr") {  ## isometric log-ratio transformation
            y <- log(w) - mean(log(w))
            y <- as.vector(y %*% t(helm(3)))
        }
        a <- numeric(nu)
        for (l in 1:nu) {
          a[l] <- as.vector(t(z[l, ] - c(y[1], y[2])) %*%
          ts %*% ( z[l, ] - c(y[1], y[2])) )
        }
        can <- 1/(2 * pi) * (1/con) * sum(exp(-0.5 * a))/nu
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
        w3 <- (2 * x2[j])/sqrt(3)
        w2 <- x1[i] - x2[j]/sqrt(3)
        w1 <- 1 - w2 - w3
        w <- c(w1, w2, w3)
        if (type == "alr") y <- log(w[-3]/w[3])  ## alr transformation
        if (type == "ilr") {  ## isometric log-ratio transformation
          y <- log(w) - mean(log(w))
          y <- as.vector(y %*% t(helm(3)))
        }
        a <- numeric(nu)
        for (l in 1:nu) a[l] <- as.vector(t(z[l, ] - c(y[1], y[2])) %*% ts %*%
          (z[l, ] - c(y[1], y[2])))
        can <- 1/(2 * pi) * (1/con) * sum(exp(-0.5 * a))/nu
        if (abs(can) < Inf) {
            mat[i, j] <- can
        } else  mat[i, j] <- NA
      }
    }
  }

  contour(x1, x2, mat, col = 3)
  proj <- matrix(c(0, 1, 1/2, 0, 0, sqrt(3)/2), ncol = 2)
  da <- x %*% proj
  points(da[, 1], da[, 2])
  b1 <- c(1/2, 0, 1, 1/2)
  b2 <- c(sqrt(3)/2, 0, 0, sqrt(3)/2)
  b <- cbind(b1, b2)
  points(b[, 1], b[, 2], type = "l", xlab = " ", ylab = " ")
}

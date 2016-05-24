################################
#### Kernel density estimation of directional data with a von Mises kernel
#### Tuning the bandwidth
#### Tsagris Michail 8/2015
#### mtsagris@yahoo.gr
#### Tuning the bandwidth
################################

vmfkde.tune <- function(x, low = 0.1, up = 1) {
  ## x is the data
  x <- as.matrix(x)  ## makes sure x is a matrix
  x <- x/sqrt(rowSums(x^2))  ## makes sure x is directional data
  p <- ncol(x)  ## dimensionality of the data
  n <- nrow(x)  ## sample size of the data
  d <- crossprod(t(x))
  diag(d) <- NA  ##  we do not want to take the diagonal elements
   funa <- function(h) {
    A <- d/h^2
    cpk <- ( (1/h^2)^(p/2 - 1) )/( (2 * pi)^(p/2) * besselI(1/h^2, p/2 - 1) )
    f <- rowSums( exp(A + log(cpk)), na.rm = TRUE )/(n - 1)
    mean(log(f))
  }
  a <- optimize(funa, c(low, up), maximum = TRUE)
  res <- c(a$maximum, a$objective)
  names(res) <- c("Optimal h", "cv")
  res
}

########################
#### Multivariate kernel density estimation
#### Tsagris Michail 2/2015
#### mtsagris@yahoo.gr
#### References: Arsalane Chouaib Guidoum (2015)
#### Kernel Estimator and Bandwidth Selection for Density
#### and its Derivatives. The kedd package
#### http://cran.r-project.org/web/packages/kedd/vignettes/kedd.pdf
#### M.P. Wand and M.C. Jones (1995)
#### Kernel smoothing, pages 91-92.
#### B.W. Silverman (1986)
#### Density estimation for statistics and data analysis, pages 76-78.
################################

mkde <- function(x, h, thumb = "none") {
  x <- as.matrix(x)  ## makes sure x is a matrix
  ## h is the h you want, which is either a vector or a single number
  ## thumb can be either "none" so the specified h is used, or
  ## "scott", or "silverman"
  n <- nrow(x)
  d <- ncol(x)  ## sample and dimensionality of x

  if (thumb == "silverman") {
    s <- apply(x, 2, sd)
    h <- (4/(d + 2))^(1/(d + 4)) * s * n^(-1/(d + 4))

  } else  if (thumb == "scott") {
    s <- apply(x, 2, sd)
    h <- s * n^(-1/(d + 4))
  } else if (thumb == "estim") {
    h <- mkde.tune(x)$hopt
  }

  else  h <- h

  if (length(h) == 1) {
    h <- diag(h, d)
  } else h <- diag(h)

  con <- prod( diag(h) )
  y <- x %*% solve(h)
  a1 <- dist(y, diag = TRUE, upper = TRUE)
  a1 <- as.matrix(a1)
  f <- as.vector( 1/(2 * pi)^(d/2) * (1/con) * rowMeans( exp(-0.5 * a1^2 ) ) )
  f
}

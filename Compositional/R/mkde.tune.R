########################
#### Multivariate kernel density estimation: tuning the bandwidth h
#### using the maximum likelihood cross validation
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

mkde.tune <- function( x, low = 0.1, up = 3, s = cov(x) ) {
  ## x contains the multivariate data
  ## low and up are the limits within which the
  ## search is conducted
  x <- as.matrix(x)
  n <- nrow(x)
  d <- ncol(x)  ## sample and dimensionality of x
  s <- s ## can put a robust covariance matrix here if you want
  ds <- 1/det(s)^0.5
  a2a <- matrix( nrow = n, ncol = n )
   tune <- function(h) {
    for (i in 1:n)  a2a[i, ] <- mahalanobis(x, x[i,], s)
    a <- a2a/h^2
    f <- ds / (2 * pi)^(d/2) * (1/h^d) *
    ( rowSums( exp(-0.5 * a) ) - 1 ) / (n - 1)
    mean( log(f) )
   }
  low <- low    ;   up <- up
  bar <- optimize(tune, c(low, up), maximum = TRUE)
  list( hopt = bar$maximum, maximum = bar$objective )
}

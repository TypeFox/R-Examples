################################
#### alpha-transformation
#### Tsagris Michail 5/2013
#### References: Tsagris M. T., Preston, S. and Wood A. T. A. (2011).
#### A data-based power transformation for
#### compositional data. In Proceedings of the 4rth Compositional Data Analysis Workshop, Girona, Spain.
#### mtsagris@yahoo.gr
################################

alfa <- function(x, a, h = TRUE) {
  ## x contains the compositional data
  ## a is the power parameter, usually between -1 and 1
  ## if h is TRUE the multiplication with the Helmert matrix takes place
  x <- as.matrix(x)  ## makes sure x is a matrix
  D <- ncol(x) ## number of components
  if ( D == 1 )   x <- t(x)
  x <- x / rowSums(x)  ## makes sure x is compositional data
  if (a != 0) {
    z <- ( x^a )
    ta <- rowSums(z)
    z <- z / ta
    z <- (D/a) * z - 1/a
    sa <- sum( log(ta) )
  } else {  ## if a=0 the ilr is calculated
    xa <- log(x)
    z <- xa - rowMeans( xa )  ## this is the clr
    sa <- nrow(x) * log(D)
  }
  if (h == TRUE) {
    aff <- z %*% t( helm( D ) ) ## multiply by the Helmert sub-matrix
    res <- list(sa = sa, aff = aff)
  } else {
    res <- list(sa = sa, aff = z)
  }
  res
}

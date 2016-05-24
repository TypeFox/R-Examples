################################
#### Frechet mean
#### Tsagris Michail 5/2013
#### References: Tsagris, M. T., Preston, S., and Wood, A. T. A. (2011).
#### A data-based power transformation for
#### compositional data. In Proceedings of the 4rth Compositional Data Analysis Workshop, Girona, Spain.
#### mtsagris@yahoo.gr
################################

frechet <- function(x, a) {
  ## x contains the compositional data
  ## a is the power parameter, usually between -1 and 1
  x <- as.matrix(x)  ## makes sure x is a matrix
  x <- x / rowSums(x)  ## makes sure x is compositional data
  if (a == 0) {
     xa <- log(x)
     y <- xa - rowMeans(xa)
     m1 <- colMeans(y)
     m <- exp(m1) / sum( exp(m1) )  ## closed geometric mean
  }  else {
     xa <- x^a
     z <- xa / rowSums(xa)
     m1 <- colMeans(z) ^ ( 1 / a )
     m <- m1 / sum(m1)  ## frechet mean in general
  }
  m
}

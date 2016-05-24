################################
### Inverse of the alpha transformation
#### Tsagris Michail 5/2013
#### References: Tsagris, M. T., Preston, S., and Wood, A. T. A. (2015).
#### Improved classification for compositional data using the alpha-transformation
#### http://xxx.tau.ac.il/pdf/1506.04976v2.pdf
#### mtsagris@yahoo.gr
################################

alfainv <- function(x, a, h = TRUE) {
  ## x is the data, not compositional
  ## a is the power parameter
  x <- as.matrix(x)
  D <- ncol(x)
  if ( D == 1)   x <- t(x)
  if (h == TRUE)  {
    h <- helm( D + 1 )  ## multiply with the Helmert
    ## sub-matrix to bring them onto Q^D
    y <- x %*% h
  }	 else y <- x
  if (a != 0) {
    z <- ( a * y + 1 )^( 1/a )
    z <- z / rowSums(z)
  } else {
    ## is a=0, the inverse of the clr is calculated
    ey <- exp(y)
    z <- ey / rowSums( ey )
  }
  z
}

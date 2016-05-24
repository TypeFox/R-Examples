################################
#### alpha-distance
#### Tsagris Michail 5/2013
#### References: Tsagris, M. T., Preston, S., and Wood, A. T. A. (2011).
#### A data-based power transformation for
#### compositional data. In Proceedings of the 4rth Compositional Data Analysis Workshop, Girona, Spain.
#### mtsagris@yahoo.gr
################################

alfadist <- function(x, a) {
  ## x contains the compositional data
  ## a is the power parameter, usually between -1 and 1
  x <- as.matrix(x)  ## makes sure x is a matric
  y <- alfa(x, a, h = TRUE)$aff
  disa <- fields::rdist(y)
  disa
}

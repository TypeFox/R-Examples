################################
#### Ridge regression with compositional data
#### in the covariates side using the alpha-transformation
#### Plot to see how things go
#### Tsagris Michail 2/2016
#### mtsagris@yahoo.gr
#### References:Tsagris M. (2015)
#### Regression analysis with compositional data containing zero values
#### Chilean journal of statistics 6(2): 47-57
################################

alfaridge.plot <- function(y, x, a, lambda = seq(0, 5, by = 0.1) ){
  ## y is dependent univariate variable. It can be matrix or vector
  ## x are compositional data, the covariates
  ## a is the value of a for the alpha-transformation
  ## lambda contains a grid of values of the ridge regularization parameter
  z <- alfa(x, a, h = TRUE)$aff ## apply the alpha-transformation
  ridge.plot(y, z, lambda = seq(0, 5, by = 0.1) )
}

################################
#### Multivariate or univariate regression with compositional data
#### in the covariates side using the alpha-transformation
#### Tsagris Michail 8/2015
#### mtsagris@yahoo.gr
#### References:Tsagris M. (2015)
#### Regression analysis with compositional data containing zero values
#### Chilean journal of statistics 6(2): 47-57
################################

alfa.pcr <- function(y, x, a, k, xnew = NULL) {
  ## y is dependent univariate variable. It can be matrix or vector
  ## x are compositional data, the covariates
  ## a is the value of a for the alpha-transformation
  ## k is the number of principal components to use
  ## oiko can be either "normal", "binomial" or "poisson"
  ## depending on the type of the independent variable
  ## "normal" is set by default
  z <- alfa(x, a, h = TRUE)$aff ## apply the alpha-transformation

  if ( length(unique(y)) == 2 ) {
    oiko <- "binomial"
  } else if ( sum(y) - round(y) == 0 ) {
    oiko <- "poisson"
  } else oiko <- "normal"

  if (oiko == 'normal') {
    mod <- pcr(y, z, k, xnew = xnew)
  } else  mod <- glm.pcr(y, z, k, xnew = xnew)
  mod ## principal component regression with the alpha-transformed
  ## compositional data
}

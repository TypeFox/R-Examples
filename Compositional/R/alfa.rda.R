################################
#### Regularised discriminant analysis for compositional data using the alpha-transformation
#### Tsagris Michail 7/2015
#### References: Tsagris, M., Preston S. and Wood A.T.A. (2016).
#### Improved classication for compositional data using the alpha-transformation
#### Journal of Classification (To appear)
#### http://arxiv.org/pdf/1506.04976v2.pdf
#### mtsagris@yahoo.gr
################################

alfa.rda <- function(xnew, x, ina, a, gam = 1, del = 0) {
  ## xnew is the new compositional observation
  ## x contains the compositional data
  ## ina is the grouping variable
  ## a is the value of the alpha
  ## gam is between pooled covariance and diagonal
  ## gam*Spooled+(1-gam)*diagonal
  ## del is between QDA and LDA
  ## del*QDa+(1-del)*LDA
  y <- alfa(x, a)$aff  ## apply the alpha-transformation
  ynew <- alfa(xnew, a)$aff
  rda(xnew = ynew, x = y, ina = ina, gam = gam, del = del)
}

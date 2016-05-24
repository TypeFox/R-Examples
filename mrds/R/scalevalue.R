# Compute scale of detection function
#
# Uses a log link
#  key.scale scale parameters
#  z design matrix for scale covariates
#
# returns Vector of scale values
# documented in ?distpdf
scalevalue <- function(key.scale, z){
  exp(as.matrix(z) %*% key.scale)
}

###############################################
###The function "FisherZ" does fisher's Z-transformation on an estimator of X, X takes values in 0~1.
###after Fisher's Z transformation, Z is normally distributed
#############################################

FisherZ <-
function(x)
  {
    
    x <- as.numeric(x)
    Z <- 0.5*log((1+x)/(1-x))
    Z
  }


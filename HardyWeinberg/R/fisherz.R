fisherz <- function(r) {
  # Fisher's z transformation for a correlation coefficient.
   y <- 0.5*log((1+r)/(1-r))
   return(y)
}

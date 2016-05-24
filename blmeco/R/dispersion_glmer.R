# Function to test for overdispersion in any model
#
# source: https://stat.ethz.ch/pipermail/r-sig-mixed-models/2011q1/015392.html
dispersion_glmer <- function(modelglmer){
  # computing  estimated scale  ( binomial model) following  D. Bates :
  # That quantity is the square root of the penalized residual sum of
  # squares divided by n, the number of observations, evaluated as:
  n <- length(resid(modelglmer))
  return(  sqrt( sum(c(resid(modelglmer),modelglmer@u) ^2) / n ) ) 
} 
#should be between, 0.75 and 1.4 if not under- or overdispersed, respectively

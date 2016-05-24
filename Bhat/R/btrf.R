"btrf" <-
function(xt,xl,xu) {

####  back transformation
####  this assumes logit transformations of the parameters
####  bounded from below by xl and from above by xu

  rho <- exp(xt)
  return((rho * xu + xl)/(1.+rho))
}

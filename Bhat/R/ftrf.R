"ftrf" <-
function(x,xl,xu) {

####  forward transformation
####  this assumes logit transformations of the parameters
####  bounded from below by xl and from above by xu
  if(any((x-xl) <= 0)) stop('ftrf requires x > xl')
  if(any((xu-x) <= 0)) stop('ftrf requires x < xu')
  return(log((x-xl)/(xu-x)))
}

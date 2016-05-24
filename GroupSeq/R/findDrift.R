"findDrift" <-
function(n,t,t2,lowerBounds,upperBounds,confidenceLevel,drift, nMax)
{
  ######################
  #INITIALIZE VARIABLES#
  ######################
  drift<-0

  ## Starting value for drift: there may be better choices.
  drift[1] <- ( upperBounds[n]+qnorm(confidenceLevel) ) / sqrt(t[n])
  drift <- computeDrift(n,t,t2,lowerBounds,upperBounds,confidenceLevel,drift[1],nMax)

  return(drift)

}#end <--*function(...)*

"computeConfidenceIntervall" <-
function(confidenceLevel, Zvalue, n, t, t2, lowerBounds, upperBounds,nMax)
{
  ###INITIALIZE VARIABLES###
  confidenceLimit<-0
  zcrit<-0

  ##Save upperBounds[n]
  tempUpperBounds <- upperBounds[n]
  upperBounds[n] <- Zvalue

  ##Use naive limits as starting values
  zcrit <- qnorm(1-(1-confidenceLevel)/2)
  confidenceLimit[1] <- (upperBounds[n]-zcrit)/sqrt(t[n])
  confidenceLimit[2] <- (upperBounds[n]+zcrit)/sqrt(t[n])

  for(i in 1:2)
  {
    ##compute target Zvalue which function 'computeDrift' has to fullfill
    target <- (i-1)*confidenceLevel + (1-confidenceLevel)/2
    confidenceLimit[i] <- computeDrift(n,t,t2,lowerBounds,upperBounds,target,confidenceLimit[i],nMax)
  }

  ##Re-Replace upperBounds[n]
  upperBounds[n] <- tempUpperBounds

  return(confidenceLimit)


}#end <--*function(...)*

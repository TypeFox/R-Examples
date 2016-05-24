"stdDeviations" <-
function( n, t2 )
{
  ###Initialize variables###
  standardDeviation<-0
  standardDeviationProcess<-0

  ##Loop for standard deviations###

  #Standard deviation for first analysis:
  standardDeviation[1]<-sqrt(t2[1])
  standardDeviationProcess[1]<-standardDeviation[1]


  #catch n<2
  if (n<2)
  {
    #do nothing else
    #user had chosen just one interim analysis
  }
  else
  {
    #Standard deviation for second and later analyses:
    for (i in 2:n)
    {
      standardDeviation[i] <- sqrt(t2[i]-t2[i-1])
      standardDeviationProcess[i] <- sqrt(t2[i])
    }
  }

  #prepare return arguments
  toBeReturned <- list(standardDeviation=standardDeviation, stdvOfTheProcess=standardDeviationProcess)
  return(toBeReturned)

}#end <--*function(...)*

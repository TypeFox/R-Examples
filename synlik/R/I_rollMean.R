# Very simple and stupid rolling mean
.rollMean <- function(vett, Lag)
{
  nObs <- length(vett)
  output <- numeric(nObs)
  
  for(ii in 1:nObs)
    output[ii] <- mean(vett[max(1, ii - Lag + 1):ii])
  
  return( output )
}



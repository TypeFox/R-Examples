percentile <- function(x, distn, include.observed=FALSE)
{
  if(include.observed)
    distn <- c(x, distn)
  
  n <- length(distn)
    
  return(findInterval(x, distn[order(distn)]) / n)
}

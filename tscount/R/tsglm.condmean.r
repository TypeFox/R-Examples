tsglm.condmean <- function(link=c("identity", "log"), ...){
  #Recursion for the conditional mean and its derivatives of a count time series following GLMs

  ##############
  #Checks and preparations:
  link <- match.arg(link)
  
  if(link == "identity") result <-  ingarch.condmean(...)
  if(link == "log") result <-  loglin.condmean(...)
  
  return(result)
}
 
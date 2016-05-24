tsglm.loglik <- function(link=c("identity", "log"), ...){
  #Conditional log-likelihood function, score function and information matrix of a count time series following GLMs

  ##############
  #Checks and preparations:
  link <- match.arg(link)
  
  if(link == "identity") result <-  ingarch.loglik(...)
  if(link == "log") result <-  loglin.loglik(...)
  
  return(result)
}
 
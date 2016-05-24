summary.bab <- function(object, ...){
  rval <- object
  class(rval) <- c("babSummary")
  return(rval)
}

summary.cab <- function(object, ...){
  rval <- object
  class(rval) <- c("cabSummary")
  return(rval)  
}


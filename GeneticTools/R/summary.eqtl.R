# Version: 30-11-2012, Daniel Fischer

summary.eqtl <- function(object, ...){
  
  ws <- object$windowSize
  if(is.null(ws)) ws <- "trans-eQTL" 

  cat("EQTL Summary\n")
  cat("---------------\n")
  cat("Type of test        :",object$method,"\n")
  cat("Tested genes        :",length(object$eqtl),"\n")
  cat("Window size (in MB) :",ws,"\n")
  invisible(object)
} 

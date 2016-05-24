##' @export
print.IPCW <- function(x,digits=3,...){
  cat("\nEstimated inverse of the probability of censoring weights (IPCW)\n\n")
  method=switch(x$method,
    "cox"="Cox regression",
    "marginal"="Kaplan-Meier",
    "nonpar"="Stratified Kaplan-Meier",
    "aalen"="Additive Aalen regression",
    "none"="No weighting",
    "Dont know")
  cat("Method for estimation: ", method,"\n")
  cat("Handler function: ",paste(as.character(x$fit$call[1]),"()",sep=""),"\n")
  if (!is.null(x$IPCW.times)){
    cat("\nhead() of the predicted IPCW for", NROW(x$IPCW.times),"subjects (rows), at the",NCOL(x$IPCW.times),"requested times (columns):\n\n") 
    print(head(x$IPCW.times),digits=digits,quote=FALSE)
  }
  if (!is.null(x$IPCW.subjectTimes)){
    cat("\nhead() of predicted IPCW at the individual subject times:\n\n")
    print(head(x$IPCW.subjectTimes),digits=digits,quote=FALSE)
  }
}

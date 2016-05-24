summary.MLoef <- function(object, ...)
{
# summary method for object of class "MLoef" (MLoef)

  # prepare message for split criteria
  if( length(object$splitcr) == 1){
    if( (object$splitcr == "median") | (object$splitcr == "mean") ){
      spl <- object$splitcr
    }
  } else {
    spl <- "user-defined"
  }

#  if(!is.null(object$warning)){
#    if(object$splitcr == "median") cat("Warning: Item(s)",paste(names(object$warning),collapse=", "),"with raw score equal to the median assigned to the lower raw score group!\n")
#    if(object$splitcr == "mean") cat("Warning: Item(s)",paste(names(object$warning),collapse=", "),"with raw score equal to the mean assigned to the lower raw score group!\n")
#  }

  cat("\n")
  cat("Martin-Loef-Test (split criterion: ",spl,")\n",sep="")

  for(i in 1:length(object$i.groups)){
    cat("\n")
    cat("Group ",i,":\nItems: ",sep="")
    cat(paste(object$i.groups[[i]]),sep=", ")
    cat("\nLog-Likelihood:",round(object$subModels[[i]]$loglik,3),"\n")
  }

  cat("\n")
  cat("Overall Rasch-Model:\n")
  cat("Log-Likelihood:",round(object$fullModel$loglik,3),"\n")
  cat("\n")

  cat(paste("LR-value:",round(object$LR,3),"\n"))
  cat(paste("Chi-square df:",round(object$df,3),"\n"))
  cat(paste("p-value:",round(object$p.value,3)),"\n")
  cat("\n")

}

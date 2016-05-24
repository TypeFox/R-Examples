print.MLoef <- function(x,...)
{
#print method for object of class "MLoef" (MLoef)

# prepare message for split criteria
  if( length(x$splitcr) == 1){
    if( (x$splitcr == "median") | (x$splitcr == "mean")){ spl <- x$splitcr }
  }
  else{ spl <- "user-defined" }
#
#  if(!is.null(x$warning)){
#    if(x$splitcr == "median") cat("Warning: Item(s)",paste(names(x$warning),collapse=", "),"with raw score equal to the median assigned to the lower raw score group!\n")
#    if(x$splitcr == "mean") cat("Warning: Item(s)",paste(names(x$warning),collapse=", "),"with raw score equal to the mean assigned to the lower raw score group!\n")
#  }
  cat("\n")
  cat("Martin-Loef-Test (split criterion: ",spl,")\n",sep="")

  cat(paste("LR-value:",round(x$LR,3),"\n"))
  cat(paste("Chi-square df:",round(x$df,3),"\n"))
  cat(paste("p-value:",round(x$p.value,3)),"\n")
  if (!("MLx" %in% class(x))) cat("\n") # no blank line if called from print.MLobj (NPtest)
}

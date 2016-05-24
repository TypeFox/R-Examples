summary.blastRes <- function(object, ...){
  cat("BlastRes Summary\n")
  cat("----------------\n")
  cat("Total Blast Results :",length(object$RID),"\n")
  cat("Used Database       :",object$database,"\n")
  invisible(object)
} 

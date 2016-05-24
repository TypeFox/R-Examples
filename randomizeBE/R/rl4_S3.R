###############################################################################
# S3 methods for the class rl4 (randomisation list)
# 
# Author: dlabes
###############################################################################

print.rl4 <- function(x, sumry=FALSE, ...){
  cat ("Randomization table          created: ")
  cat(format(x$date),"\n")
  cat ("(seed:", x$seed, "blocksize:", x$blocksize,")\n\n")
  print(x$rl, row.names = FALSE, ...)
  # summary if wished
  if(sumry) {
    cat("\n\nSummary of randomisation\n\n")
    cat(nrow(x$rl), "subjects randomized into ")
    cat(nrow(unique(x$rl["sequence"])), "sequence groups.\n")
    cat("Number of subjects in sequence groups:\n")
    print(x$ninseqs)
    pval <- x$runs.pvalue
    if(!is.null(pval)) {
      cat("Runs test of randomness: p.value=")
      if(pval>0.0001){
        cat(formatC(pval, format="f", digits=4),"\n")
      } else {
        cat(" <0.0001\n")
      }
    }      
  }
  return(invisible(x))
}


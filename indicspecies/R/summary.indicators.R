summary.indicators <- function (object, ...) {
  cat(paste("Result of 'indicators()':\n\n"))
  cat(paste("  Number of plots in target site group: ", sum(object$group.vec),"\n",sep=""))
  cat(paste("  Number of candidate species: ", length(object$candidates),"\n",sep=""))
	cat(paste("  Number of species used in combinations: ", length(object$candidates),"\n",sep=""))
  cat(paste("  Number of indicators (single species or species combinations) kept: ", nrow(object$C),"\n",sep=""))
  cat(paste("  Group coverage: ", coverage(object)*100,"%\n",sep=""))
  cat("\n")
}

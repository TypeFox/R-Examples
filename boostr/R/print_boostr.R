#' @export
print.boostr <- function(x, ...) {
  reweighterOutput <- reweighterOutput(x)
  numEstimators <- length(ensembleEstimators(x))
  performance <- estimatorPerformance(x)
  
  cat(paste0("A boostr object composed of ", numEstimators, " estimators.\n"))
  cat(paste0("\nAvailable performance metrics: ",
             paste0(names(performance[[1]]), collapse=", "), "\n"))
  cat("\nStructure of reweighter output:\n")
  cat(str(reweighterOutput), "\n")
  if (extractCalcBoostrPerformance(x)) {
    cat("Performance of Boostr object on Learning set:\n")
    print(extractPerformanceOnLearningSet(x))
    cat("\n")
  }
  invisible(x)
}
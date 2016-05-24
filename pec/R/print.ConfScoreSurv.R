##' @export
print.confScoreSurv <- function(x,...){
  overall <- do.call("cbind",lapply(x$models,function(m){
    colMeans(m$score)
  }))
  mm <- cbind(x$times,overall)
  colnames(mm) <- c("times",names(x$models))
  rownames(mm) <- rep("",NROW(mm))
  cat("\nPopulation average confidence score:\n\n")
  print(mm)
}

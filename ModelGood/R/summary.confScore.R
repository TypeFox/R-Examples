summary.confScore <- function(object,...){
  print(object$splitMethod)
  cat("\n")
  avscore <- do.call("rbind",lapply(object$models,function(x){
    mean(x$score)
  }))
  colnames(avscore) <- "MeanConfScore"
  print(avscore)
}

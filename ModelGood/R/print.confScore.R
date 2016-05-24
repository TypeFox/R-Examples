print.confScore <- function(x,...){
  print(x$splitMethod)
  cat("\n")
  avscore <- do.call("rbind",lapply(x$models,function(x){
    mean(x$score)
  }))
  colnames(avscore) <- "MeanConfScore"
  print(avscore)
}

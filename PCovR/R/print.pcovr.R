print.pcovr <-
function(x, ...){
  cat("Loading matrix:\n")
  srt <- SortLoadings(x$Px)
  print(t(srt))
  cat("\nRegression weight matrix:\n")
  print(t(x$Py))
  cat("\nExplained variance in X:",x$Rx2,"\n")
  cat("\nExplained variance in Y:",x$Ry2,"\n")
  if (x$modsel=="sim"){
    cat("\nCross-validation fit:\n")
    print(x$Qy2)
  }
  if (x$modsel=="seq"){
    cat("\nWeighted variance accounted for:\n")
    print(x$VAFsum)
  }
  if (x$modsel=="seqAcv"){
    cat("\nWeighted variance accounted for:\n")
    print(x$VAFsum)
    cat("\nCross-validation fit:\n")
    print(x$Qy2)
  }
  if (x$modsel=="seqRcv"){
    cat("\nCross-validation fit:\n")
    print(x$Qy2)
  }
}

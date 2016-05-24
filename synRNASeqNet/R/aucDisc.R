aucDisc <-
function(fpr, tpr){
  f <- cbind(fpr, tpr)
  f <- f[order(f[, 1], f[, 2]), ] 
  ans <- 0
  fpr <- f[, 1]
  tpr <- f[, 2]
  ans <- sapply(2:length(fpr), function(i){
    (fpr[i] - fpr[i-1])*tpr[i-1] + 0.5*(fpr[i] - fpr[i-1]) * (tpr[i] - tpr[i-1])})
  ans <- sum(ans)
  return(ans)
}

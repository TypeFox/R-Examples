Profile.ICC <-
function(set1, set2, ..., omit=TRUE) {

  mat <- abind(set1, set2, ..., along=3)
  if(omit==T) {
    res <- t(apply(mat, 1, function(x) get.ICC(na.omit(x))))
  }
  if(omit==F) {
    res <- t(apply(mat, 1, get.ICC))
  }
  colnames(res) <- c("ICC1", "ICC1k", "ICC2", "ICC2k", "ICC3", "ICC3k")
  return(res)
}

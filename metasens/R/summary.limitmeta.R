summary.limitmeta <- function(object, ...){

  meta:::chkclass(object, "limitmeta")
  
  res <- object
  class(res) <- c("summary.limitmeta", "limitmeta")
  
  res
}

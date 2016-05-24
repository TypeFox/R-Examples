composite <- function(set, R = NULL, Zitems = FALSE, maxScore = NULL, rel=FALSE, nomiss = 0.8, tr = 0) {
  if(Zitems) {
    set <- scale(set)
    if(!is.null(R)) {
      revSet <- set[,R]*-1
      set <- data.frame(set[,-R], revSet)
    }
  }
  if(!Zitems) {
    if(!is.null(R)) {
      if(is.null(maxScore)) {
        maxScore <- max(set)
      }
      revSet <- set[,R]
      revScored <- 1 + maxScore - revSet
      set <- data.frame(set[,-R], revScored)
    }
  }
  if(rel) {
    rel.stats <- alpha(set)$total
    print(rel.stats)
  }
  return(sapply(data.frame(t(set)), meanif, nomiss, tr))
}
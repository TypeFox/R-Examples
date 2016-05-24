CompPVs <-
function(T){
  tmp <- sort(T)
  ind <- order(T)
  m <- length(T)
  pv <- rep(1, length(T))
  i <-1
  for(j in seq_len(m-1)) {
    if (tmp[j] < tmp[j + 1]) {
      pv[ind[i:j]] <- j / m
      i <- j + 1
    }		
  }
  pv
}


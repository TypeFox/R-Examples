freq.abs <- function(data) {
  #Function produces absolute frequencies count for every variable in a dataset.
  r <- nrow(data)
  s <- ncol(data)
  
  freq.table <- matrix(data = 0, ncol = s, nrow = max(unique(data)))
  
  for (i in 1:s) {
    for (j in 1:max(data[ ,i])) {
      count <- vector(mode="numeric", length=0)
      for (k in 1:r) {
        if (data[k,i] == j) {
          count[k] <- 1
        }
        else {
          count[k] <- 0
        }
      }
      freq.table[j,i] <- sum(count)
    }
  }
  return(freq.table)
}
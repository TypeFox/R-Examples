initWeights <- function(diss) {
  if (!is.list(diss)) {
	  n <- attr(diss,"Size")
    ww <- matrix(1, n, n)
    ww[is.na(as.matrix(diss))] <- 0                ## blank out missings
	  return(as.dist(ww))
  } else {
  n <- attr(diss[[1]],"Size")
  m <- length(diss)
  ww <- repList(matrix(1,n,n),m)
  for (i in 1:m) {
    wwi <- ww[[i]]
    wwi[is.na(as.matrix(diss[[i]]))] <- 0
    ww[[i]] <- as.dist(wwi)
  }
  return(ww)
  }
}


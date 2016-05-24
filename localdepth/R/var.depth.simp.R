#############################################################
#
#	var.depth.simp function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: June, 23, 2008
#	Version: 0.1
#
#	Copyright (C) 2008 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

var.depth.simp <- function(x, nsamp='exact') {
  x <- as.matrix(x)
  nr <- nrow(x)
  nc <- ncol(x)
  nt <- choose(nr, nc+1)
  if (is.numeric(nsamp) && nsamp <= 0) stop("the argument 'nsamp' must be positive")
  if (is.numeric(nsamp) && nsamp > nt) {
      warning("Since 'nsamp' is greater than the number of simplex the 'exact' method is used")
      nsamp <- 'exact'
  }

  if (is.character(nsamp) && nsamp=='exact') {
    nsamp <- nt
    res <- rep(0, nsamp)
    pos <- (nc+1):1
    div <- nr-(0:nc)
    for (i in 1:(nsamp-1)) {
      res[i] <- volume.simp(x[pos,])
###      cat(i, pos, '\n')
      temp <- pos%%div
      temp0 <- rev((1:(nc+1))[temp==0])
      tempno0 <- (1:(nc+1))[temp!=0]
      pos[min(tempno0)] <- pos[min(tempno0)] + 1
      if (length(temp0)) {
        for (j in 1:length(temp0)) {
          pos[temp0[j]] <- max(pos[(temp0[j]+1):(nc+1)])+1
        }
      }
    }
    res[nsamp] <- volume.simp(x[pos,])
  } else if (is.numeric(nsamp)){
    res <- rep(0, nsamp)
    for (i in 1:nsamp) {
      index <- sample((1:nr), size=(nc+1), replace=FALSE)
      res[i] <- volume.simp(x[index,])
    }
  } else {
    stop("the argument 'nsamp' must be either 'exact' or a positive number")
  }
  result <- list()
  result$depth.var <- mean(res)
  result$volumes <- res
  return(result)
}


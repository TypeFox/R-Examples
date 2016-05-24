### nPairs.R
###------------------------------------------------------------------------
### What: Number of variable pairs - code
### $Id$
### Time-stamp: <2008-12-30 18:29:58 ggorjan>
###------------------------------------------------------------------------

nPairs <- function(x, margin=FALSE, names=TRUE, abbrev=TRUE, ...)
{
  ## --- Setup ---
  if(!is.data.frame(x) & !is.matrix(x)) stop("'x' must be a data.frame or a matrix")
  k <- ncol(x)
  if(!margin) {
    ret <- matrix(nrow=k, ncol=k)
  } else {
    ret <- matrix(nrow=k, ncol=k + 1)
  }

  ## --- Count ---
  diag(ret)[1:k] <- apply(X=x, MARGIN=2, FUN=function(x) sum(!is.na(x)))
  for(i in 1:k) {
    for(j in i:k) {
      ret[i, j] <- ret[j, i] <- sum(!is.na(x[, i]) & !is.na(x[, j]))
      if(margin) {
        if(i == 1) {
          ret[i, (k + 1)] <- ret[1, 1]
        } else {
          ret[i, (k + 1)] <- sum(rowSums(!is.na(x[, c(1:i)])) == i)
        }
      }
    }
  }
  
  ## --- Names ---
  if(names) {
    tmp <- colnames(x)
    if(abbrev) tmp <- as.character(abbreviate(tmp, ...))
    rownames(ret) <- tmp
    if(margin) {
      colnames(ret) <- c(tmp, "all")
    } else {
      colnames(ret) <- tmp
    }
  }
  class(ret) <- c("nPairs", class(ret))
  ret
}

summary.nPairs <- function(object, ...)
{
  n <- nrow(object)
  ret <- matrix(data=0, nrow=n, ncol=n)
  for(i in 1:n) {
    tmp <- 1:n
    tmp <- tmp[!(tmp == i)]
    ret[i, tmp] <- object[i, i] - object[i, tmp]
  }
  dimnames(ret) <- dimnames(object)
  ret
}

###------------------------------------------------------------------------
### nPairs.R ends here

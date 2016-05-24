## This file (R/truthTable.R) is part of QCA3 package
## copyright: HUANG Ronggui 2008-2012

coincid <- function(x,standardized=FALSE,use=c("complete","pairwise")){

  coincid_fn <- function(x,y,standardized,...){
    ## helper function
    allvalues <- !is.na(x) & !is.na(y)
    x<-x[allvalues]
    y<-y[allvalues]
    Min <- pmin(x,y)
    Max <- pmax(x,y)
    Sumx <- sum(x)
    Sumy <- sum(y)
    if (standardized){
      ans <- sum(Min)/min(Sumx,Sumy)
    } else {
      ans <- sum(Min)/sum(Max)
    }
    return(ans)
  }

  minmax <- range(x,na.rm=TRUE)
  if (minmax[1] < 0 || minmax[2] > 1) stop("All the values of 'x' should be range from 0 to 1.")
  use <- match.arg(use)
  if (use=="complete") x= na.exclude(x)
  Nvar <- ncol(x)
  ans <- matrix(numeric(0),nrow=Nvar,ncol=Nvar)
  index <- t(combn(Nvar,2)) # the first col is the column index
  nindex <- nrow(index)
  for (i in 1:nindex) {
    rindex <- index[i,][2]
    cindex <- index[i,][1]
    ans[rindex,cindex] <- coincid_fn(x[,rindex],x[,cindex],standardized=standardized)
  }
  rownames(ans) <- colnames(ans) <- names(x)
  diag(ans) <- 1
  class(ans) <- "coincid"
  ans
}

print.coincid <- function(x,digits=3,...)
{
  cat("Coincidence Matrix\n\n")
  x<-unclass(x)
  print(x,digits=digits,na.print=" ",quote = FALSE,...)
}

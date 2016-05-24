which.maxn <- function(x,n=1){
  if (n==1)
    which.max(x)
  else
    {
      if (n>1){
        ii <- order(x,decreasing=TRUE)[1:min(n,length(x))]
        ii[!is.na(x[ii])]
      }
      else {
       stop("n must be >=1")
      }
    }
}

which.minn <- function(x,n=1){
  if (n==1)
    which.min(x)
  else
    {
      if (n>1){
        ii <- order(x,decreasing=FALSE)[1:min(n,length(x))]
        ii[!is.na(x[ii])]
      }
      else {
       stop("n must be >=1")
      }
    }
}



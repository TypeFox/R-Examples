mlag <- function (x, lag) {
  len=NROW(x)
  ncol=NCOL(x)
  index=matrix(rep(0:(len-1), ncol), nrow=len) 
  lag=lag %% len
  index=(index-matrix(rep(lag, each=len), nrow=len)) %% len + 1
  result=matrix(nrow=len, ncol=ncol, NA)
  for (i in 1:ncol)
    result[, i]=x[index[, i], i]  
  return(result)
}

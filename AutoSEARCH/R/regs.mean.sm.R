regs.mean.sm <-
function(y, mc=NULL, ar=NULL, ewma=NULL, mx=NULL)
{
y.n <- length(y)
mX <- NULL

#constant
if(!is.null(mc)){mX <- cbind(mX,rep(1,y.n)); colnames(mX) <- "mconst"}

#ar terms:
arlags <- NULL
if(!is.null(ar)){
  for(i in 1:length(ar)){
    arlags <- cbind(arlags, gLag(y, k=ar[i]))
  }
  colnames(arlags) <- paste("ar", ar, sep="")
}
mX <- cbind(mX, arlags)

#EqWMA term:
if(is.null(ewma)){EqWMA <- NULL}else{
  EqWMA <- do.call(eqwma, c(list(y),ewma) )
}
mX <- cbind(mX, EqWMA)

#create matrix of mean regressors mx:
if(!is.null(mx)){
  mxnames <- colnames(mx)
  if(is.null(mxnames)){
    mxnames <- paste("mx", 1:ncol(mx), sep="")
  }
  if(any(mxnames == "")){
    missing.colnames <- which(mxnames == "")
    for(i in 1:length(missing.colnames)){
      mxnames[i] <- paste("mx", i, sep="")
    }
  }
  colnames(mx) <- mxnames
}
mX <- cbind(mX, mx)

#out-matrix:
out <- mX

return(out)

}

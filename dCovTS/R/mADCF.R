mADCF <- function(x,lags,output=TRUE){
 if(!is.matrix(x)) x <- as.matrix(x)
 if ((is.data.frame(x)) | (is.matrix(x))){
  if(NCOL(x)==1)stop('Only multivariate time series with dimension d>2')
  if(!is.ts(x)) x <- as.ts(x)
 }
 if (!is.numeric(x)) 
     stop("'x' must be numeric")
 if(!all(is.finite(x))) stop('Missing or infitive values')
 if (missing(lags)) 
       stop("'lags' is missing with no default")
 n <- as.integer(NROW(x))
 d <- as.integer(NCOL(x))
 if ((lags >= 0) && (lags <= (n-1))){
  X <- rbind(x[(1+lags):n,])
  Y <- rbind(x[1:(n-lags),])
 } else if ((lags >= (-n+1)) && (lags <0)){
  X <- rbind(x[1:(n+lags),])
  Y <- rbind(x[(1-lags):n,])
 } else {
  stop("lags must be in the range of -(n-1) and (n-1)")
 }
 cross.adcf <- matrix(NA,d,d)
 for (i in 1:d){
  for (j in 1:d){
   cross.adcf[i,j] <- dcor(X[,i],Y[,j])
  }
 }
 if(output){
 cat("Distance Correlation Matrix at lag: ", lags, "\n")
 print(cross.adcf)
 }
 mADCF <- cross.adcf
}


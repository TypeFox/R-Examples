##create a data matrix
t <- seq(0,4*pi,len=50)
X.org <- cbind(cos(t),sin(2*t)) %*% matrix(rnorm(20),2,10)

##add some normal errors
X <- X.org + .25*rnorm(length(X.org))
##and mark some data as missing
X[runif(length(X))<.25] <- NA

##Ensure that we have complet columns/rows
while( any(rowSums(is.na(X))==dim(X)[2]) || any(colSums(is.na(X))==dim(X)[1]) ){
  X <- X.org + .25*rnorm(length(X.org))
  X[runif(length(X))<.25] <- NA
}

#run the missing data svd
res <- SVDmiss(X, niter=100, ncomp=2)
#look at the status
res$status
#plot the first four columns of the data matrix
par(mfrow=c(2,2))
for(i in 1:4){
  plot(t, X[,i])
  lines(t, res$Xfill[,i])
  lines(t, X.org[,i], col=2)
}
\dontshow{
  if( max(abs(res$Xfill-X),na.rm=TRUE) > 1e-10 ){
    stop("SVD.miss: Points not interpolated")
  }
}

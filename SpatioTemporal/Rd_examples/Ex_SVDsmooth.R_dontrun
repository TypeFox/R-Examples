##create a data matrix
t <- seq(0,4*pi,len=50)
X.org <- cbind(cos(t),sin(2*t)) %*% matrix(rnorm(10),2,5)

##add some normal errors
X <- X.org + .25*rnorm(length(X.org))
##and mark some data as missing
X[runif(length(X))<.25] <- NA

##Ensure that we have complet columns/rows
while( any(rowSums(is.na(X))==dim(X)[2]) || any(colSums(is.na(X))==dim(X)[1]) ){
  X <- X.org + .25*rnorm(length(X.org))
  X[runif(length(X))<.25] <- NA
}

##compute two smooth basis functions
res <- SVDsmooth(X, n.basis=2, niter=100)

##or compute the function that gives the basis functions
res.fnc <- SVDsmooth(X, n.basis=2, niter=100, fnc=TRUE)

##and they are equal
summary( res.fnc()-res )
\dontshow{
  if( !isTRUE(all.equal(res.fnc(), res)) ){
    stop("SVDsmooth: Function and matrix NOT equal.")
  }
}

##plot the two smooth basis functions
par(mfcol=c(3,2), mar=c(4,4,.5,.5))
plot(t, res[,1], ylim=range(res), type="l")
lines(t, res[,2], col=2)
##and some of the data fitted to the smooths
for(i in 1:5){
  plot(t, X[,i])
  lines(t, predict.lm(lm(X[,i]~res), data.frame(res)) )
  lines(t, X.org[,i], col=2)
}

##compute cross-validation for 1 to 4 basis functions
res.cv <- SVDsmoothCV(X, n.basis=0:4, niter=100)

##study cross-validation results
print(res.cv)
summary(res.cv)

##plot cross-validation statistics
plot(res.cv, sd=TRUE)
##boxplot of CV statistics for each column
boxplot(res.cv)
##plot the BIC for each column
plot(res.cv, "BIC", pairs=TRUE)

\dontshow{
  res.cv.2 <- SVDsmoothCV(X, n.basis=0:4, niter=100, fnc=TRUE)
  if( !isTRUE(all.equal(res.cv[1:6], res.cv.2[1:6])) ){
    stop("SVDsmoothCV: Function and matrix NOT equal; stats!")
  }
  
  tmp <- lapply(res.cv.2$smoothSVD[[3]], function(x){x()})
  if( any(sapply(1:length(tmp),
                 function(i){max(abs(res.cv$smoothSVD[[3]][,,i] - tmp[[i]]))})
          >1e-14) ){
    stop("SVDsmoothCV: Function and matrix NOT equal; parts!")
  }
}

localshape <- function(xdata,proportion=0.1,mscatter="mcd",mcdalpha=0.8,
                       covstandard="det"){
#  if (mscatter=="mcd") require(robustbase)
  xdata <- as.matrix(xdata)
  scatter <- switch(mscatter,
                    mcd=covMcd(xdata,alpha=mcdalpha)$cov,
                    cov=cov(xdata))
  n <- nrow(xdata)
  p <- ncol(xdata)
  np <- round(proportion*n)
  mmatrix <- matrix(0,n,n)
  for (i in 1:n)
    mmatrix[i,] <- mahalanobis(xdata,xdata[i,],scatter)
  lcov <- matrix(0,p,p)
  for (i in 1:n){
    xc <- cov(xdata[order(mmatrix[i,])[1:np],])
    lcov <- lcov+switch(covstandard,
                        trace=xc/sum(diag(xc)),
                        det=xc/det(xc),
                        none=xc)    
  }
  lcov <- lcov/n
  lcov
}
    

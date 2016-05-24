
irf<-function(data, bParams, item, plotERF=TRUE, 
              thetaEAP = NULL, 
              minCut = -3, maxCut = 3, NCuts = 9){
## irf can be used to plot model-implied and empirical item response functions for item sets that fit
##          different polynomial IRT models


  if(is.data.frame(bParams)) bParams <- as.matrix(bParams)
  
  if(ncol(bParams)!=9) stop("\n\nbParams should have 9 columns")
 
  x <- seq(-4,4,by=.01)
  Numx <- length(x)
  xpoly <- matrix(cbind(1,x, x^2, x^3, x^4, x^5, x^6, x^7),Numx, 8)

  # Prob of a keyed response 2PL
  P <- function(m){
    1/(1+exp(-m))
  }
  
# plot k=1 ICC
plot(x,P(xpoly %*% bParams[item,1:8]), 
     typ="l",
     lwd=2,
     xlim=c(-4,4), ylim=c(0,1),
     main=paste("Item ",item, "            k = ", bParams[item,9], sep="" ),
     xlab=expression(theta), cex.lab=1.3,
     ylab="Probability")

if(!is.null(thetaEAP)){
   if(plotERF){
     erfOUT<- erf(thetaEAP, data, whichItem=item, min=minCut, max=maxCut,Ncuts=NCuts)
     points(erfOUT$centers, erfOUT$probs, pch=16, col="red", cex=1)
   }
}

}#END irf
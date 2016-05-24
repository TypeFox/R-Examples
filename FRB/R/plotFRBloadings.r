plotFRBloadings <- function(x, confmethod = c("BCA","basic"), pcs=1:min(5, length(x$eigval)), nvars=min(10, length(x$eigval))) {

FRBres <- x
confmethod <- match.arg(confmethod)
conf <- FRBres$conf
currentAsk <- devAskNewPage(ask = NULL)

q <- length(FRBres$eigval)
if (any(!(pcs %in% c(1:q)))) stop(paste("indices in 'pcs' should be in 1:",q,sep=""))

# pass the variable names to the eigenvectors
if (!is.null(dimnames(FRBres$est$Mu)[[2]]))
    dimnames(FRBres$eigvec) <- list(dimnames(FRBres$est$Mu)[[2]], paste("PC",1:q,sep=""))
else
    dimnames(FRBres$eigvec) <- list(paste("V",1:q,sep=""), paste("PC",1:q,sep=""))
    
par(mfrow=c(1,1))

for (comp in pcs) {
  PChere <- FRBres$eigvec[,comp]
  if (confmethod=="basic") {
    PChere.low <- FRBres$eigvec.CI.basic[((comp-1)*q+1):(comp*q),1]
    PChere.high <- FRBres$eigvec.CI.basic[((comp-1)*q+1):(comp*q),2]
  }
  else {
    PChere.low <- FRBres$eigvec.CI.bca[((comp-1)*q+1):(comp*q),1]
    PChere.high <- FRBres$eigvec.CI.bca[((comp-1)*q+1):(comp*q),2]
  }
  
  orderPC <- order(abs(PChere), decreasing=TRUE)
  varsTP <- PChere[orderPC[1:nvars]]
  namesvarsTP = names(PChere)[orderPC[1:nvars]]
  varsTPhigh <- PChere.high[orderPC[1:nvars]]
  varsTPlow <- PChere.low[orderPC[1:nvars]]
  
  plot(1:nvars, varsTP, type="n", ylim=c(-1.0,1.0), xaxt="n", xlab="", xlim=c(1-0.25, nvars+0.25), cex.axis=1.5, ylab="loadings", main=paste("Loadings PC",comp, " (+ ",conf*100, "% ",confmethod," confidence limits)",sep=""))
  axis(side=1, at = 1:nvars, labels = namesvarsTP)
#  grid(ny=NA, nx=nvars, lty=1)
  for (k in 1:nvars) { lines(c(k,k), c(-1,1), col="grey") }
  grid(nx=NA, ny=NULL, lty=2)
  points(1:nvars, varsTP, pch=20, col="red", cex=2)
  abline(h=-1, lwd=2)
  abline(h=1, lwd=2)
  abline(h=0, lwd=2)
  #points(1:nvars, varsTPhigh, pch=20, cex=2)
  #points(1:nvars, varsTPlow, pch=20, cex=2)
  for (i in 1:nvars) {
    lines(c(i-0.2, i+0.2), c(varsTPlow[i],varsTPlow[i]), lwd=2)
    lines(c(i-0.2, i+0.2), c(varsTPhigh[i],varsTPhigh[i]), lwd=2)
    lines(c(i-0.2, i-0.2),c(varsTPlow[i], varsTPlow[i]+0.02), lwd=2)
    lines(c(i+0.2, i+0.2),c(varsTPlow[i], varsTPlow[i]+0.02), lwd=2)
    lines(c(i-0.2, i-0.2),c(varsTPhigh[i], varsTPhigh[i]-0.02), lwd=2)
    lines(c(i+0.2, i+0.2),c(varsTPhigh[i], varsTPhigh[i]-0.02), lwd=2)
  }
  devAskNewPage(ask = TRUE) 
}
devAskNewPage(ask = currentAsk)

}

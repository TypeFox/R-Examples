HWGenotypePlot <- function(X,plottype=1,xlab=expression(f[AA]),ylab=ifelse(plottype==1,expression(f[AB]),
    expression(f[BB])),asp=1,pch=19,xlim=c(0,1),ylim=c(0,1),cex=1,cex.axis=2,cex.lab=2,...) {
  # Makes a scatter plot of genotype frequencies.
if (is.vector(X)) {
   if (length(X) != 3) {
      stop("X must have three elements")
   }
   else {
      X <- matrix(X, ncol = 3, dimnames = list(c("1"), 
      names(X)))
   }
}
nr <- nrow(X)
nc <- ncol(X)
if (any(X < 0)) 
    stop("X must be non-negative")
if (nc != 3) 
    stop("X must have three columns")
if (nrow(X) == 1) {
    Xcom <- X/sum(X)
   } else {
      Xcom <- HWClo(X)
   }
   fAA <- seq(0,1,by=0.01)
   fAB <- 2*(sqrt(fAA)-fAA)
   fBB <- (1-sqrt(fAA))^2   
   if(is.element(plottype,c(1,2))) {
     opar <- par(mar=c(5,5,2,1))
     if(plottype==1) { # heterozygote versus homozygote
       plot(Xcom[,1],Xcom[,2],pch=pch,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,cex=cex,cex.axis=cex.axis,cex.lab=cex.lab,...)
       lines(c(0,1),c(1,0),lwd=2,col="red")     
       points(fAA,fAB,pch=19,col="blue",type="l",lwd=2)
     }
     if(plottype==2) { # homozygote versus homozygote
       plot(Xcom[,1],Xcom[,3],pch=pch,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,cex=cex,cex.axis=cex.axis,cex.lab=cex.lab,...)
       lines(c(0,1),c(1,0),lwd=2,col="red")     
       points(fAA,fBB,pch=19,col="blue",type="l",lwd=2)
     }
     par(opar)
   } else stop("HWGenotypePlot: invalid argument for plottype")
return(NULL)
}

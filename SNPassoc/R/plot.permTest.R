plot.permTest<-function(x,...)
 {

  if (!inherits(x, "permTest")) 
      stop("x must be an object of class 'permTest'")

  param<-x$param
  pmin<-x$pmin
  psig<-x$psig

  o<-density(pmin) 
  grid<-seq(0,max(o$x),length=1000)
  plot(grid,dbeta(grid,param[1],param[2]),xlab="minimum p value",ylab="density",type="n")
  hist(pmin,prob=TRUE,col="gray90",border="gray50",add=TRUE)
  lines(grid,dbeta(grid,param[1],param[2]),col="blue",lwd=2)

  segments(psig,0,psig,dbeta(psig,param[1],param[2]),col="red")
  legend("topright",c("empirical distribution","theoretical distribution",paste("adjusted p value:",round(psig,8))),
        lty=c(1,1,1),col=c("gray90","blue","red"),title="",cex=0.8,bty="n")

 }

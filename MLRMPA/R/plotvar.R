plotvar <-
function(descriptor.tr,activity.tr,maintitle,xlab,breaks){

  if(missing(maintitle)){
    maintitle="Histogram of the Variance of the Descriptors"
  }
  if(missing(xlab)){
    xlab="Variance of the Descriptors"
  }
  if(missing(breaks)){
    breaks=20
  }
  
  op<-par(cex.lab=1.2)
  x<-apply(descriptor.tr,2,var)
  y<-var(activity.tr)
  tmp01<-hist(x,breaks=breaks,plot=F)
  tmp02<-tmp01$counts
  tmp02[which(tmp02==0)]<-NA
  
  hist(x,xlab=xlab,main=maintitle,breaks=breaks,ylim=c(0,1.2*max(tmp01$counts)))
  text(tmp01$mids,tmp02,labels=tmp02,cex=0.8,font=2,pos=3)
  box()
  points(y,0)
  par(op)
}

plotcor <-
function(descriptor.tr,activity.tr,maintitle,xlab,breaks,xlim){
 
  if(missing(maintitle)){
    maintitle="Histogram of the Correlation between Descriptors and Activity"
  }
  if(missing(xlab)){
    xlab="Correlation between Descriptors and Activity"
  }
  if(missing(breaks)){
    breaks=20
  }
  x<-cor(descriptor.tr,as.numeric(activity.tr[,1]))
  if(missing(xlim)){
    tmp<-pmax(abs(min(x)),max(x))
    xlim=c(-tmp,tmp)
  }
  
  op<-par(cex.lab=1.2)
  tmp01<-hist(x,breaks=breaks,plot=F)
  tmp02<-tmp01$counts
  tmp02[which(tmp02==0)]<-NA

  hist(x,xlab=xlab,
       main=maintitle,breaks=breaks,xlim=xlim,xaxt="n",ylim=c(0,1.2*max(tmp01$counts)))
  axis(1,xaxp=c(-round(tmp,1),round(tmp,1),round(tmp,1)/0.05))
  text(tmp01$mids,tmp02,labels=tmp02,cex=0.8,font=2,pos=3)
  box()
  par(op)
}

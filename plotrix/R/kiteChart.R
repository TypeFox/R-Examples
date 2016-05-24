kiteChart<-function(x,xlim=NA,ylim=NA,timex=TRUE,main="Kite chart",
 xlab=ifelse(timex,"Time","Groups"),ylab=ifelse(timex,"Groups","Time"),
 border=par("fg"),col=NULL,varpos=NA,varlabels=NA,varscale=FALSE,
 timepos=NA,timelabels=NA,mar=c(5,4,4,4),axlab=c(1,2,3,4),
 normalize=FALSE,shownorm=TRUE,...) {

 # leave a bit more space on top if there is an axis
 if(varscale || normalize) mar<-mar+c(0,0,!timex,timex)
 dimx<-dim(x)
 if(normalize) {
  if(is.na(varpos[1])) varpos<-1:dimx[2]
  kitewidths<-rep(1,dimx[2])
  kitemax<-1
 }
 else {
  if(is.na(varpos[1])) {
   kitewidths<-apply(as.matrix(x),1,max,na.rm=TRUE)
   varpos<-rep(0,length(kitewidths)-1)
   varpos[1]<-1.1*kitewidths[1]/2
   for(kite in 2:length(kitewidths))
    varpos[kite]<-varpos[kite-1]+1.1*(kitewidths[kite-1]+kitewidths[kite])/2
   kitemax<-1.1*sum(kitewidths)
  }
  else {
   kitemax<-max(diff(varpos))
   kitewidths<-apply(as.matrix(x),1,max,na.rm=TRUE)
  }
 }
 oldmar<-par(mar=mar)
 if(is.na(xlim[1])) {
  if(timex) xlim<-c(1,dimx[2])
  else {
   if(normalize) xlim<-c(0.5,dimx[1]+0.5)
   else xlim<-c(0,kitemax)
  }
 }
 if(is.na(ylim[1])) {
  if(timex) {
   if(normalize) ylim<-c(0.5,dimx[1]+0.5)
   else ylim<-c(0,kitemax)
  }
  else ylim<-c(1,dimx[2])
 }
 plot(0,xlim=xlim,ylim=ylim,main=main,xlab=xlab,ylab=ylab,type="n",
  axes=FALSE,...)
 if(is.na(varpos[1])) varpos<-1:dimx[1]
 if(is.na(varlabels[1])) {
  if(is.null(rownames(x))) varlabels<-varpos[1:dimx[1]]
  else varlabels<-rownames(x)
 }
 axis(ifelse(timex,axlab[2],axlab[1]),at=varpos[1:dimx[1]],labels=varlabels)
 if(is.na(timepos[1])) timepos<-1:dimx[2]
 if(is.na(timelabels[1])) {
  if(is.null(colnames(x))) timelabels<-timepos
  else timelabels<-colnames(x)
 }
 axis(ifelse(timex,axlab[1],axlab[2]),at=timepos,labels=timelabels)
 if(varscale && !normalize) {
  plotlim<-par("usr")
  mtext(round(kitewidths,1),side=3+timex,at=varpos)
  axis(3+timex,at=c(varpos-kitewidths/2,varpos+kitewidths/2),
   labels=rep("",2*length(varpos)))
 }
 box()
 if(is.null(col)) col<-rainbow(dimx[1])
 if(length(col) < dimx[1]) col<-rep(col,length.out=dimx[1])
 for(krow in 1:dimx[1]) {
  if(normalize) {
   if(shownorm)
    mtext(paste("*",signif(1/max(x[krow,]),digits=3)),
     ifelse(timex,axlab[4],axlab[3]),at=varpos[krow],las=1)
   x[krow,]<-x[krow,]/(max(x[krow,]))
  }
  xpos<-1:length(x[krow,])
  if(timex)
   polygon(c(xpos,rev(xpos)),
    c(varpos[krow]+x[krow,]/2,
    varpos[krow]-rev(x[krow,])/2),
    col=col[krow],border=border)
  else
   polygon(c(varpos[krow]+x[krow,]/2,
    varpos[krow]-rev(x[krow,])/2),
    c(xpos,rev(xpos)),
    col=col[krow],border=border)
 }
 invisible(oldmar)
}

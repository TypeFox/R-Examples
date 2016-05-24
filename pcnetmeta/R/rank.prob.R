rank.prob <- function(nma.obj,cex.axis=1,cex.lab=1){
  if(is.null(nma.obj$TrtRankProb)) stop("users did not specify rank.prob in the argument param of the function which produces nma.obj.")
  if(!is.null(nma.obj$TrtRankProb)){
    rank.prob<-nma.obj$TrtRankProb
  }
  ntrt<-dim(rank.prob)[1]
  trtname<-rownames(rank.prob)
  par(mar=c(4,4,0.5,1.5)+0.1,mgp=c(2.5,0.5,0))
  plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab="Treatments",ylab="Rank Probabilities",cex.lab=cex.lab)
  axis(2,at=seq(0,1,0.1),labels=TRUE,cex.axis=cex.axis)
  axis(1,at=(9/20/ntrt+seq(from=0,by=1/ntrt,length.out=ntrt))*(1-4/(10*ntrt)),labels=trtname,tick=FALSE,cex.axis=cex.axis)
  axis(4,at=c(1/(2*ntrt),1-1/(2*ntrt)),labels=c("No. 1",paste("No.",ntrt)),tick=FALSE,cex.axis=cex.axis,pos=1)
  rgb.val<-seq(from=0,by=0.8/(ntrt-1),length.out=ntrt)
  for(i in 1:ntrt){
    trt.i.rank.prob<-rank.prob[i,]
    cum.prob<-c(0,cumsum(trt.i.rank.prob))
    xleft.i<-((i-1)/ntrt)*(1-4/(10*ntrt))
    xright.i<-(i/ntrt-1/(10*ntrt))*(1-4/(10*ntrt))
    for(j in 1:ntrt){
      rect(xleft=xleft.i,ybottom=cum.prob[j],xright=xright.i,ytop=cum.prob[j+1],col=rgb(rgb.val[j],rgb.val[j],rgb.val[j]),border="black")
    }
  }
  for(i in 1:ntrt){
    rect(xleft=1-2/(10*ntrt),ybottom=(i-1)/ntrt,xright=1,ytop=i/ntrt,col=rgb(rgb.val[i],rgb.val[i],rgb.val[i]),border="black")
  }
}
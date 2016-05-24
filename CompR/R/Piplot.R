Piplot<-function(Pi,SigmaPi=NULL,level=0.05,main=NULL,ylab="Bradley's scores",xlab="Item",labelprod=NULL)
{
  names<-factor(rownames(Pi),levels=rownames(Pi),ordered=FALSE)
  if (level<=0 | level>=1)
  {
    SigmaPi<-NULL
    print("No confidence Interval, incorrect value for the level")
  }
  
  if (is.null(SigmaPi))
  {
    plot(Pi,pch="-",ylim=c(-0.05,max(1/(length(Pi)),Pi+0.1)),cex=1.5,ylab=ylab,xlab=xlab,xaxp=c(0,nrow(Pi),nrow(Pi)),xaxt="n",main=main)
  }
  else
  {
    plot(Pi,pch="-",ylim=c(-0.05,max(1/(length(Pi)),max(Pi[,1]+1.96*SigmaPi)+0.1)),cex=1.5,ylab=ylab,xlab=xlab,xaxp=c(0,nrow(Pi),nrow(Pi)),xaxt="n",main=main)
    for (i in 1:nrow(Pi))
    {
      IC<-matrix(nrow=2,ncol=2)
      IC[,1]<-i
      a<-qnorm(1-0.5*level, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
      IC[1,2]<-Pi[i]-a*SigmaPi[i]
      IC[2,2]<-Pi[i]+a*SigmaPi[i]
      arrows(IC[1,1],IC[1,2],IC[2,1],IC[2,2],lwd=1,angle=89.9,code=3,length=0.05)
    }
  }
if (is.null(rownames(Pi)))
{
  if (is.null(labelprod))
  {
    numprod<-as.character(c(1:length(Pi)))
    prefix<-rep("P",length(Pi))
    labelprod<-paste(prefix,numprod,sep="")
  }
}
else
{
  labelprod<-rownames(Pi)
}
  text(c(1:nrow(Pi)),-0.05,labels=labelprod)
}
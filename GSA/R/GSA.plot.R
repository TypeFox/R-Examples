GSA.plot=function(GSA.obj,  fac=1, FDRcut=1){

eps=.001
a=GSA.listsets(GSA.obj, geneset.names=NULL,  FDRcut=FDRcut)
neg=matrix(as.numeric(as.character(a$negative[,-2])),ncol=4)+eps
pos=matrix(as.numeric(as.character(a$positive[,-2])),ncol=4)+eps

ymax=max(c(neg[,4],pos[,4]))
xmax=max(c(neg[,3],pos[,3]))
ymin=min(c(neg[,4],pos[,4]))
xmin=min(c(neg[,3],pos[,3]))

o1=order(neg[,3])
o2=order(pos[,3])
plot(jitter(neg[o1,3],factor=fac),jitter(neg[o1,4],factor=fac),xlab="p-value",ylab="False discovery rate", type="n",log="xy", xlim=c(xmin,xmax),
ylim=c(ymin,ymax))
points(jitter(neg[o1,3],factor=fac),jitter(neg[o1,4],factor=fac),col=3,cex=.7,type="b")
points(jitter(pos[o2,3],factor=fac),jitter(pos[o2,4],factor=fac),col=4, cex=.7,type="b")
#axis(3, at = res[,3], lab = paste(round(res[,3],0)), srt = 90, adj = 0,cex=.5)
legend(.1,.2,c("Negative","Positive"),col=3:4, pch=c("o","o"))
}



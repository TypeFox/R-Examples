plotparti<-function(pa,d1=NULL,d2=NULL,
dendat=NULL,restri=NULL,pch=21,support=pa$support,col="black",cex.axis=1)
{
# support=NULL

recs<-pa$recs
if (!is.null(d1)) recs<-recs[,c(2*d1-1,2*d1,2*d2-1,2*d2)]

xmin<-min(recs[,1])
xmax<-max(recs[,2])
ymin<-min(recs[,3])
ymax<-max(recs[,4])

if (!is.null(dendat)) plot(dendat,xlab="",ylab="",pch=pch,cex.axis=cex.axis)
else plot(0,0,type="n",ylim=c(ymin,ymax),xlab="",ylab="",xlim=c(xmin,xmax),
     cex.axis=cex.axis)

len<-dim(recs)[1]

if (is.null(restri)) restric<-matrix(1,len,1)
else{
  restric<-matrix(0,len,1)
  restrilen<-length(restri)
  for (i in 1:restrilen){
     restric[restri[i]]<-1
  }
}

i<-1
while (i<=len){

 if (restric[i]==1){
    lines(c(recs[i,1],recs[i,1]),c(recs[i,3],recs[i,4]),col=col)
    lines(c(recs[i,1],recs[i,2]),c(recs[i,4],recs[i,4]),col=col)
    lines(c(recs[i,2],recs[i,2]),c(recs[i,4],recs[i,3]),col=col)
    lines(c(recs[i,2],recs[i,1]),c(recs[i,3],recs[i,3]),col=col)
 }

 i<-i+1
}

if (!is.null(support)){
  lines(c(support[1],support[1]),c(support[3],support[4]),col=col)
  lines(c(support[1],support[2]),c(support[4],support[4]),col=col)
  lines(c(support[2],support[2]),c(support[4],support[3]),col=col)
  lines(c(support[2],support[1]),c(support[3],support[3]),col=col)
}

}








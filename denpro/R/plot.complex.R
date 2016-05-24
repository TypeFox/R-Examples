plot.complex<-function(complex,dendat,xlab="",ylab="",cex.lab=1,cex.axis=1,pch=19,
col=NULL,border="black")
{
plot(dendat,xlab=xlab,ylab=ylab,cex.lab=cex.lab,cex.axis=cex.axis,pch=pch)
lkm<-dim(complex)[1]
for (i in 1:lkm){
    cur<-complex[i,]
    x<-dendat[cur,1]
    y<-dendat[cur,2]
    polygon(x,y,col=col,border=border)
}

}


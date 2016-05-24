plot.HiddenF <-
function(x,y=NULL,main="Hidden Additivity Plot",rfactor="Rows Factor",cfactor="Columns Factor",colorvec=c("black","red"),legendx=FALSE,center=FALSE,...)
{
hfobj <- x
ymtx <- makemtx.fcn(hfobj$tall)
a <- nrow(ymtx)
b <- ncol(ymtx)
grpvector <- 1+hfobj$config.vector[1+b*(0:(a-1))]
#matplot(t(ymtx),type="l",col=grpvector,main=main,xaxt="n",ylab="y",xlab=cfactor,lty=grpvector,...) # xlab=cfactor
if(center==TRUE)
{
ymtx <- ymtx-apply(ymtx,1,mean)
}
matplot(t(ymtx),type="l",col=colorvec[grpvector],main=main,xaxt="n",ylab="y",xlab=cfactor,lty=1:a,...) # xlab=cfactor
if (legendx == TRUE){
legend(x=locator(1),legend=1:a,lty=1:a,col=colorvec[grpvector],title=rfactor)}
axis(1,at=1:b,labels=names(table(hfobj$tall$cols)),cex.axis=1.2)
}

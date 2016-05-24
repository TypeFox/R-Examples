`ogive.freq` <-
function(histogram,type="",xlab="",ylab="",axes="",las=1,...)
{
if(axes=="") ejes=TRUE
else {
ejes<-FALSE
if (axes) ejes<-TRUE
}
if(type=="") type<-"b"
#if(las=="") las<-1
    yy <- histogram$counts
    zz <- histogram$breaks
    y1 <- sum(yy)
    nx <- length(yy)
    zz <-c(zz,2*zz[nx+1]-zz[nx])
    nz<-length(zz)
    y<-rep(0,nz)
    for (i in 1:nx) {
        y[i+1] <- y[i] + yy[i]/y1
    }
    y[nz]<-1
probability<-y
if(ejes){
plot(zz,y,type=type,xlab=xlab,ylab=ylab,axes=FALSE,...)
axis(1,zz,las=las)->ax; axis(2,las=las)->ay
abline(v=ax,h=ay,lty=2,col="grey")
}
else plot(zz,y,type=type,xlab=xlab,ylab=ylab,axes=axes,...)
table<-data.frame(x=zz,RCF=round(y,4))
if(xlab=="")xlab<-"x"
names(table)[1]<-xlab
return(table)
}



paracoor<-function(X,Y=NULL,xmargin=0.1,
paletti=matrix("black",dim(X)[1],1),noadd=TRUE,verti=NULL,cex.axis=1,
points=TRUE,col.verti="black",col.verti.y="red",digits=3,
arg=NULL,colarg="red",lwd=1,cex=1,yaxt="s")
{
n<-dim(X)[1]
d<-dim(X)[2]
ylim<-c(min(X),max(X))
if (is.null(Y)) D<-d else D<-d+dim(Y)[2]

if (noadd)
plot(x="",y="",
xlim=c(1-xmargin,D+xmargin),ylim=ylim,
xlab="",ylab="",xaxt="n",cex.axis=cex.axis,yaxt=yaxt)

for (i in 1:n){
    if (points) points(X[i,],col=paletti[i],cex=cex)
    for (j in 1:(d-1)) segments(j,X[i,j],j+1,X[i,j+1],
                                col=paletti[i],lwd=lwd)
}
#if (points) for (i in 1:n) points(X[i,],col=paletti[i])

if (!is.null(Y)){
  miny<-min(Y)
  maxy<-max(Y)
  z<-matrix(0,n,dim(Y)[2])
  for (i in 1:n){
     for (j in 1:dim(Y)[2]){
         coeff<-(Y[i,j]-miny)/(maxy-miny)
         z[i,j]<-ylim[1]+coeff*(ylim[2]-ylim[1])
     }
  }
  for (i in 1:n){
     j<-2
     while (j<=dim(Y)[2]){ 
        if (points) points(d+j,z[i,j],col=paletti[i],cex=cex)  
        segments(d+j-1,z[i,j-1],d+j,z[i,j],col=paletti[i],lwd=lwd)
        j<-j+1
     }
     if (points){
         points(d+1,z[i,1],col=paletti[i],cex=cex)
         points(d,X[i,d],col=paletti[i],cex=cex)
     }
     segments(d,X[i,d],d+1,z[i,1],col=paletti[i],lwd=lwd)
  }
  #if (points) for (i in 1:n) points(d+1,z[i],col=paletti[i])
  segments(d+0.5,ylim[1],d+0.5,ylim[2],col=col.verti.y,lwd=lwd)
  text(d+dim(Y)[2]+xmargin/2,ylim[1],format(miny,digits=digits))
  text(d+dim(Y)[2]+xmargin/2,ylim[2],format(maxy,digits=digits))
  text(d+dim(Y)[2]+xmargin/2,ylim[1]+(ylim[2]-ylim[1])/2,
       format(miny+(maxy-miny)/2,digits=digits))
}

if (!is.null(verti)) segments(verti,ylim[1],verti,ylim[2],col=col.verti,lwd=lwd)

if (!is.null(arg)){
    if (points) points(arg,col=colarg,cex=cex)
    for (j in 1:(d-1)) segments(j,arg[j],j+1,arg[j+1],col=colarg,lwd=lwd)
}

}


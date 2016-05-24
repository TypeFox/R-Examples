kmeans.assig.groups=function(out,draw=TRUE,...){
if (!is.null(out$lcenters))  lxm=out$lcenters
else  lxm=NULL
mdist=out$z.dist
nr=nrow(out$fdataobj)
nc=ncol(out$fdataobj)
xm=out$centers[["data"]]
ncl=nrow(xm)
grupo=rep(0,nr)
d=out$d
par(mfrow=c(1,2))
ngroups=nrow(d)-nrow(out$fdataobj[["data"]])
for (i in 1:nr){    grupo[i]=which.min(d[(nr+1):(nr+ngroups),i])      }
if (draw){
 if (nr==2){
  plot(out$fdataobj,main="Assigning groups")
  for (i in 1:ngroups){points(xm[i,1],xm[i,2],col=i+1,pch=8,cex=1.5)}
  }
 else{
  plot(out$fdataobj,col=grupo+1,main="Assigning groups",lwd=.3,lty=2)
  lines(out$centers,col=2:(length(grupo+1)),lwd=3,lty=1)     #new
 }
}
if (nc==2) { for (j in 1:nc){points(xm[j,1],xm[j,2],col=j+1,pch=7,cex=1.5)}}
return(list("centers"=out$centers,"cluster"=grupo))
}
##







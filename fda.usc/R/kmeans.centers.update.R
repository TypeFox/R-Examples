
kmeans.centers.update=function(out,group,dfunc=func.trim.FM,draw=TRUE,par.dfunc=list(trim=0.05),...){
if (class(out)!="kmeans.fd") stop("Error: incorrect input data")
z=out$fdataobj[["data"]]
tt=out$fdataobj[["argvals"]]
rtt<-out$fdataobj[["rangeval"]]
names=out$fdataobj[["names"]]
mdist=out$z.dist
xm=out$centers[["data"]]
centers=out$centers
nr=nrow(z)
nc=ncol(z)
grupo=group
ngroups=length(unique(group))
d=out$d
ncl=nrow(xm)
for (j in 1:ngroups){
     if (sum((grupo==j))>0) {
               dm=z[grupo==j,]
               ind=which(grupo==j)
               if (is.vector(dm) || nrow(dm)<3) {k=j}#revisar pq  k no hace nada!!
               else   {
                     par.dfunc$fdataobj<-centers
                     par.dfunc$fdataobj$data<-dm

                     stat=do.call(dfunc,par.dfunc)
                     }

               if (is.fdata(stat)) xm[j,]=stat[["data"]]
               else  xm[j,]=stat
                    }
 }
centers$data=xm
rownames(centers$data)<-paste("center ",1:ngroups,sep="")
if (draw){
 if (nr==2){
  plot(out$fdataobj,main="Center update")
  for (i in 1:ngroups){points(xm[i,1],xm[i,2],col=i+1,pch=8,cex=1.5)}}
 else{
plot(out$fdataobj,col="grey",lty=grupo+1,lwd=0.15,cex=0.2,main="Update centers")
lines(centers,col=2:(length(grupo+1)),lwd=3,lty=1)
   }}
return(list("centers"=centers,"cluster"=grupo))
}


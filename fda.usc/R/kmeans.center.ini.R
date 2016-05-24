kmeans.center.ini=function(fdataobj,ncl=2,metric=metric.lp,draw=TRUE,method="sample",iter=100,par.metric=NULL,...){
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
if (is.null(par.metric)) par.metric=list()
par.metric$fdata1<-fdataobj
#mdist=metric(fdataobj,...)
mdist=do.call(metric,par.metric)
z<-fdataobj[["data"]]
tt<-fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
names<-fdataobj[["names"]]
nr=nrow(fdataobj)
nc=ncol(fdataobj)
#else stop("The argument metric is not properly defined")
if (is.vector(ncl)) {
     len.ncl=length(ncl)
     if    (len.ncl==1) {
        ind=1
        ngroups=ncl
        if (method=="sample")  {
#            lxm=sample(nr,ngroups,replace=FALSE) #max.iter=1
              max.combn<-choose(nr,ngroups)
              max.iter<-min(iter,max.combn)
              vec<-array(NA,dim=c(ngroups,max.iter))
              vec.d<-rep(NA,max.iter)
              for (i in 1:max.iter)  {
               vec[,i]<-sample(1:nr,ngroups,replace=FALSE)
               vec.d[i]<-sum(mdist[vec[,i],vec[,i]])
               }
               ind.max<-which.max(vec.d)
               lxm<-vec[,ind.max]
              }
        else if (method=="exact") {
          co<-combn(1:nr,ngroups)
          vec<-rep(NA,ncol(co))
          for (i in 1:ncol(co)) {vec[i]<-sum(mdist[co[,i],co[,i]])}
          max.vec<-which.max(vec)
          lxm<-co[,max.vec]
          }
         else stop("Center initialization method unknown")
          xm=z[lxm,]
        }
        else      stop("the argument ncl is expected the number of groups to detect")
      }         else      stop("the argument ncl is expected the number of groups to detect")
 d=rbind(mdist,mdist[lxm,])
 centers=fdata(xm,tt,rtt,names)
if (draw){
 if (nr==2){
  plot(fdataobj)
  for (i in 1:ngroups){points(xm[i,1],xm[i,2],col=i+1,pch=8,cex=1.5)}}
 else{
  plot(fdataobj,col="grey",lty=2)
  lines(centers,col=2:(ngroups+1),lwd=3,lty=1)
#  for (i in 1:ngroups){    points(tt,xm[i,],col=i+1,lwd=3)}
  }
}
out=list("centers"=centers,"lcenters"=lxm,"z.dist"=mdist,"fdataobj"=fdataobj)
class(out)="kmeans.fd"
return(invisible(out))
}

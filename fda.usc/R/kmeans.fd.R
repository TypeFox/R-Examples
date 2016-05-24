
kmeans.fd=function(fdataobj,ncl=2,metric=metric.lp,dfunc=func.trim.FM,max.iter=100,par.metric=NULL,par.dfunc=list(trim=0.05),
par.ini=list(method="sample"),draw=TRUE,...) {
#if (is.data.frame(z)) z=as.matrix(z)
#else if (is.vector(z))     z <- as.matrix(t(z))
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
 nas1<-apply(fdataobj$data,1,count.na)
 if (any(nas1))  stop("fdataobj contain ",sum(nas1)," curves with some NA value \n")
z<-fdataobj[["data"]]
tt<-fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
nr=nrow(z)
nc=ncol(z)
if (is.vector(ncl)) {
  len.ncl=length(ncl)
  if (len.ncl==1) {
     par.ini$fdataobj=fdataobj
#      extras2 <- match.call(expand.dots = FALSE)$...
     if (is.null(par.ini$method)) par.ini$method="sample"
     if (is.null(par.ini$ncl))  par.ini$ncl=ncl
     if (is.null(par.ini$metric)) par.ini$metric=metric
     if (is.null(par.ini$draw)) par.ini$draw=draw
     if (is.null(par.ini$iter)) par.ini$iter=100
     if (!is.null(par.metric)) par.ini$par.metric<-par.metric
#      par.ini$...<-extras2
      par.ini$...<-par.metric
#     (fdataobj,ncl=2,metric=metric.lp,draw=TRUE,method="sample",iter=100,...){
     out1=do.call(kmeans.center.ini,par.ini)
#      out1=kmeans.center.ini(fdataobj=fdataobj,ncl=ncl,metric=metric,draw=draw,...)
     lxm<-out1$lcenters
     out1$d=rbind(out1$z.dist,out1$z.dist[lxm,])
     }
else {
     ngroups=length(ncl)
     lxm=ncl
     xm=z[lxm,]
     out1=list()
     out1$fdataobj<-fdataobj
     out1$ncl=len.ncl
#     mdist=metric(fdataobj,...)   #hacer do.call
if (is.null(par.metric)) par.metric=list("p"=2,"w"=1)
par.metric$fdata1<-fdataobj
#mdist=metric(fdataobj,...)
mdist=do.call(metric,par.metric)
     out1$z.dist<-mdist
     out1$d=rbind(mdist,mdist[lxm,])
     out1$centers<-fdataobj[ncl,]
     out1$lcenters<-ncl
     class(out1)="kmeans.fd"
     }}
 else if (is.fdata(ncl)) {   # fdata centers
   lxm=NULL
   xm=ncl[["data"]]
if (is.null(par.metric)) par.metric=list("p"=2,"w"=1)
par.metric$fdata1<-fdataobj
#mdist=metric(fdataobj,...)
mdist=do.call(metric,par.metric)
par.metric2<-par.metric
par.metric2$fdata2<-ncl
mdist2=do.call(metric,par.metric2)
#   mdist=metric(fdataobj,...)
#   mdist2=metric(fdataobj,ncl,...)    ##### metric
   out1=list()
   out1$fdataobj<-fdataobj
   out1$centers=ncl
   out1$lcenters<-NULL
   ngroups=nrow(ncl)
   ncl=nrow(ncl)
   out1$d=rbind(mdist,t(mdist2))
   class(out1)="kmeans.fd"
}
 ngroups=nrow(out1$centers[["data"]])
 a=0;aa<-i<-1
# C <- match.call()
# mf <- match.call(expand.dots = FALSE)
 same_centers=FALSE
#while ((i<max.iter) && (a!=aa)) {
while ((i<max.iter) && (!same_centers)) {
  out3=kmeans.assig.groups(out1,draw=draw)
  out2=kmeans.centers.update(out1,group=out3$cluster,dfunc=dfunc,draw=draw,par.dfunc=par.dfunc,...)
  a=out2$cluster
  aa=out3$cluster
  same_centers<-out2$centers$data==out3$centers$data
  out1$centers<-out2$centers
  i=i+1
     }
#cat("iterations: ",i)
out<-list("cluster"=out2$cluster,"centers"=out2$centers)
return(out)
}


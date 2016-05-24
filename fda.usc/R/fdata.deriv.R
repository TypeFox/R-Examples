fdata.deriv<-function(fdataobj,nderiv=1,method="bspline",class.out='fdata'
,nbasis=NULL,...) {
 if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
 nas1<-apply(fdataobj$data,1,count.na)
 if (any(nas1))  stop("fdataobj contain ",sum(nas1)," curves with some NA value \n")
 DATA<-fdataobj[["data"]]
 tt=fdataobj[["argvals"]]
  rtt=fdataobj[["rangeval"]]
 labels=fdataobj[["names"]]
 ndist=ncol(DATA)
 if (method=="diff") {
  res=matrix(NA,nrow=nrow(DATA),ncol=ncol(DATA))
  for (i in 1:nrow(DATA)) {
   a=diff(DATA[i,],differences=nderiv)/(tt[2:ndist]-tt[1:(ndist-1)])
   ab=matrix(NA,ncol=ndist,nrow=2)
   ab[1,2:ndist]=a
   ab[2,1:(ndist-1)]=a
#   res[i,]=apply(ab,2,mean,na.rm=TRUE)
   res[i,]=colMeans(ab,na.rm=TRUE)
  }
  labels$main<-paste("d(",labels$main,",",nderiv,")",sep="")
  res<-fdata(res,tt,rtt,names=labels)
 }
  else {
      if (any(method==c("fmm", "periodic", "natural", "monoH.FC"))) {
       res=matrix(NA,nrow=nrow(DATA),ncol=ncol(DATA))
       for (i in 1:nrow(DATA)) {
         f1<-splinefun(x=tt,y=DATA[i,],method=method)
         res[i,]=f1(tt,deriv=nderiv)
        }
          labels$main<-paste("d(",labels$main,",",nderiv,")",sep="")
         res<-fdata(res,tt,rtt,names=labels)
       }
      else{
      if (any(method==c("bspline","exponential", "fourier",
      "monomial","polynomial"))) {
#no run  "constant","polygonal","power"
#            res<-fdata(DATA,tt)
      res=fdata2fd(fdataobj=fdataobj,type.basis=method,nbasis=nbasis,nderiv=nderiv,...)
      if (class.out=='fdata') {
         ff<-eval.fd(tt,res)
         labels$ylab<-paste("d(",labels$ylab,",",nderiv,")",sep="")
         res=fdata(t(ff),tt,rtt,names=labels)
         }
      }
  }
}
res
}



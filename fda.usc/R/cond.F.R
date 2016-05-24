cond.F=function(fdata0,y0,fdataobj,y,h=0.15,g=0.15,metric=metric.lp,
Ker=list(AKer=AKer.epa,IKer=IKer.epa),...){
 if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
  nas1<-apply(fdataobj$data,1,count.na)
tt<-fdataobj$argvals
 if (any(nas1))  stop("fdataobj contain ",sum(nas1)," curves with some NA value \n")
if (!is.fdata(fdata0))  fdata0=fdata(fdata0,tt)
  nas2<-apply(fdata0$data,1,count.na)
 if (any(nas2))  stop("fdata0 contain ",sum(nas2)," curves with some NA value \n")
 if (any(is.na(y0)))  stop("y0 contain ",sum(is.na(y0)),"  NA values \n")
 if (any(is.na(y)))   stop("y contain ",sum(is.na(y)),"  NA values \n")
data<-fdataobj[["data"]]
n = nrow(data)
m = ncol(data)
nn = nrow(fdata0)
ndist=length(y0)
  xy.dist=metric(fdataobj,fdata0,...)
  x.dist=metric(fdataobj,...)
  h3=quantile(x.dist, probs = h, na.rm = TRUE,type=4)
  h=h3
  W =apply(xy.dist/h,2,Ker$AKer)
  g=g*(max(y)-min(y))
  A=outer(y,y0,"-")
  Wy = apply(A/g,2,Ker$IKer) ###
  WW= W
  if (nn==1)  {
      s=sum(WW)
      Fc =1-drop(t(Wy)%*%WW/s)
 }
 else  {
#    s=apply(WW,2,sum,na.rm=TRUE)
    s=colSums(WW,na.rm=TRUE)
    s=diag(1/s,nrow=length(s),ncol=length(s))
    Fc =1-drop((t(Wy)%*%WW)%*%s)
    }
  out=list("Fc"=Fc,"y0"=y0,"g"=g,"h"=h,"x.dist"=x.dist,"xy.dist"=xy.dist)
return(out)
}

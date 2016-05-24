sensitivity.plot<-function(y,sub,stats){
  ra<-range(y)
  xr<-mean(ra)+c(-1,1)*diff(ra)
  outlier<-xr[1]+(0:100)*diff(xr)/100
  base<-stats(y)
  sens<-array(NA,c(length(base),length(outlier)))
  dimnames(sens)[[1]]<-names(base)
  for(j in seq(length(outlier))){
    sens[,j]<-stats(c(y,outlier[j]))-base
  }
  plot(xr,range(sens),type="n",main="Sensitivity",sub=sub,
       ylab="Change in Statistic Value",xlab="New Observation")
  inds<-seq(length(base))
  for(i in inds) lines(outlier,sens[i,],lty=i,col=i)
  legend(min(ra),max(sens),lty=inds,col=inds,
         legend=names(base))
}




dmannwhitney<-function(u,m,n){
  if((u<0)|(u>(m*n))){
    val<-0
  }else{
    if(m==1) val<-1
    if(n==1) val<-1
    if((m>1)&(n>1)) val<-dmannwhitney(u,m-1,n)+dmannwhitney(u-m,m,n-1)
  }
  return(val)
}


betatest<-function(x,y){
  if(length(x)>10) cat("Warning: this will take forever.")
  out<-.Fortran("betatest"
                ,as.integer(length(x))
                ,as.double(x)
                ,as.double(y)
                ,pval=as.double(0)
                ,PACKAGE="MultNonParam"
                )
  return(out$pval)
}




wilding<-function(u1,u2,m1,n1,m2,n2){
  out<-.Fortran("wildings",u1=as.integer(u1)
                ,u2=as.integer(u2)
                ,m1=as.integer(m1)
                ,n1=as.integer(n1)
                ,m2=as.integer(m2)
                ,n2=as.integer(n2)
                ,out=as.double(0.0)
                ,PACKAGE="MultNonParam")
  return(out$out)
}


util.jplot<-function(x,y,...){
  newx<-newy<-rep(NA,length(y))
  newy[1]<-y[1]; newx[1:2]<-x[1]
  ry<-diff(range(y))
  begin<-1
  for(j in 2:length(y)){
    if(abs(y[j]-newy[begin])>(ry*.001)){ 
      begin<-begin+1; 
      newy[begin]<-y[j]; 
      newx[2*begin-(1:0)]<-x[j]
    }else{
      newx[2*begin]<-x[j]
    }
  }
  newx<-newx[seq(2*begin)]
  newy<-rep(newy[seq(begin)],rep(2,begin))
  plot(newx,newy,type="n",...)
  for(j in seq(length(newy)/2)){
    if(newx[2*j]==newx[2*j-1]){
      points(newx[2*j],newy[2*j])
    }else{
      lines(newx[2*j-(1:0)],newy[2*j-(1:0)])
    }
  }
  return(invisible(list(x=newx,y=newy)))
}


testve<-function(n,m,k,nsamp=100,delta=0,beta=0,disc=0){
  o<-rep(NA,nsamp)
  i<-rep(c(rep(0,m),rep(1,n)),k)
  str<-rep(1:k,rep(n+m,k))
  for(j in seq(nsamp)){
    x<-rnorm(k*(n+m))
    y<-rnorm(k*(n+m))
    if(disc!=0) y<-round(y/disc)*disc
    ds<-as.data.frame(list(str=str,x=x,y=y+(i-1)*delta+beta*x, i=i))
    out <-probest(ds,"y","i","str","x",0)
    o[j]<-(out$b-.5)/sqrt(out$Vb)
  }
  st<-paste("Delta=",delta,"beta=",beta,"m=",m,"n=",n)
  if(disc!=0) st<-paste(st,"made discrete")
  browser()
  aaa<-qqnorm(o,plot.it=F)
  plot(aaa$x,aaa$y,main="Normal QQ Plot for Distribution of KKW Statistic",
       sub=st,type="p")
  abline(b=1,a=0)
  #  qqline(o)
}


#ds :data set
#resp : response manifest variable vector
#grp  : vector of the variable used to form groups (Advanced/non-advanced for the prostate example, Recurrent/non-recurrent for the breast cancer example)
#str : strata variable vector
probest<-function(ds,resp,grp,str,covs=NULL,delta=NA){
  r<-length(resp)
  if(length(delta)==1){if(is.na(delta)) delta<-rep(0,r)}
  if(length(delta)==0) delta<-rep(0,r)
  covariates<-if(length(covs)<2){
    if(is.null(covs)){
      as.double(NULL) 
    }else{
      if(is.factor(ds[,covs])){
        model.matrix(~ds[,covs]-1)[,-1,drop=F]
      }else{
        ds[,covs,drop=F]
      }
    }
  }else{
    ds[,covs]
  }
  M<-as.integer(if(is.null(covs)) 0 else dim(covariates)[2])
  N<-dim(ds[,resp,drop=F])[1]
  grpv<-as.numeric(as.factor(ds[[grp]]))
  gn<-sort(unique(grpv))
  strv<-as.numeric(as.factor(ds[[str]]))
  ustr<-unique(strv)
  out<-.Fortran("probest",
                as.integer(length(resp)), M, as.integer(N), as.integer(grpv),
                as.integer(length(gn)), as.integer(gn),
                as.integer(strv), as.integer(ustr),as.integer(length(ustr)),
                as.double(as.matrix(ds[,resp])), as.double(as.matrix(covariates)),
                as.logical(!is.na(ds[,resp])),as.double(delta),
                b=as.double(rep(0,r)),Vb=as.double(rep(0,r^2))
                ,PACKAGE="MultNonParam"
  )
  return(list(b=out$b,Vb=array(out$Vb,c(r,r))))
}

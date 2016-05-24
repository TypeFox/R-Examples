
# Metro_Hastings1=function(fcn,estimates,prop_sigma,iterations=iterations){
#   mt=dim(prop_sigma)[1]
#   prop=ginv(.1*diag(mt)+ginv(prop_sigma))
#   th1=mvrnorm(1,estimates,prop)
#   trace=matrix(0,iterations,length(estimates))
#   for (i in 1:(iterations+1000)){
#     #browser()
#     th2=mvrnorm(1,th1,diag(diag(prop_sigma)))
#     ratio=exp(fcn(th2)-fcn(th1))
#     if(!is.na(ratio)) {
#       r=min(ratio,1)
#       if(rbinom(1,1,r)==1) th1=th2}
#       if(i>1000) trace[i-1000,]=th1
#     }
#   return(list(trace=trace))
#   }

pssm.noninferiority.boundary<-function(x,time,cov1,cov2,
                                       approximate=TRUE,alpha=0.05,iterations=50000){
  #   if(!is.null(delta)){
  #     if(is.null(Q)){Print('ERROR delta requires Q');return(NULL)} 
  #     else {precision=1/(-delta/qnorm((1-Q)/2))}} 
  #   else if(!is.null(Sd)){precision=1/Sd^2}
  curvf=pssm.survivalcurv(x,cov1,cov2,timeToProgression=FALSE,covariance=TRUE)
  contrast=matrix(c(1,-1),1,2)
  upper<-function(precision){
    m=length(precision)
    upper.out=rep(0,m)
    for (i in 1:m){
      curv=curvf(time,precision[i])
      variance=attributes(curv)$covariance
      sd=sqrt(contrast%*%variance%*%t(contrast))
      upper.out[i]=sum(contrast*curv$estimate)-qnorm(1-alpha)*sd}
    #browser()
    return(upper.out)
  }
  if(approximate)  return(upper)
  else {
  upper1<-function(precision){
    m=length(precision)
    upper.out=rep(0,m)
    for (i in 1:m){
    t=as.list(x@call)
    t$prior=c(0,precision)
    xx=eval(as.call(t))
    prop_sigma=xx@covariance.estimates
    mt=dim(prop_sigma)[1]
    prop=ginv(.1*diag(mt)+ginv(prop_sigma))
    mh=Metro_Hastings(xx@loglike,xx@estimates,prop,iterations=iterations)

    #mh=Metro_Hastings1(xx@loglike,xx@estimates,prop_sigma=xx@covariance.estimates,iterations=iterations)
    ms=dim(mh$trace)[1]
    st=round(ms*.25) #removing first 25% during adaption.
    value=rep(0,ms-st+1)
    j=0
    obj=x
    for (jj in st:ms){
      j=jj+1
      obj@estimates=mh$trace[j,]
      curv=pssm.survivalcurv(obj,cov1,cov2,timeToProgression=FALSE,covariance=FALSE)(time,0)
      value[j]=sum(curv$estimate*contrast)
    }
    upper.out[i]=quantile(value,alpha)
    }
    #attributes(upper.out)<-list(density=density(value),values=value)
 
    return(upper.out)

    }  
    return(upper1)
}
}

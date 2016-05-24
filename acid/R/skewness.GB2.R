skewness.GB2 <-
function(b,a,p,q){
  skewnessGB2<-rep(NA,length(b))
  for(j in 1:length(b)){
    fmGB2<- km.GB2(b[j],a[j],p[j],q[j],k=1) # first moment
    smGB2<- km.GB2(b[j],a[j],p[j],q[j],k=2) # second moment
    tmGB2<- km.GB2(b[j],a[j],p[j],q[j],k=3) # third moment
    if(-a[j]*p[j]<3&a[j]*q[j]>3){
      skewnessGB2[j] <- (tmGB2-3*fmGB2*(smGB2-fmGB2^2)-fmGB2^3)/(smGB2-fmGB2^2)^(3/2)
    } else print(paste("-ap=",-a[j]*p[j],"and aq=",a[j]*q[j]))
  }
  return(skewnessGB2)
}

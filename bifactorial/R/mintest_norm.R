#min-test for bifactorial designs (Laska and Meisner, 1989)
#Implementation taken from Hellmich and Lehmacher (2005)
pvstudent2<-function(C,nboot,simerror,...){
  if(is.binary(C@data)) warning("For binary data, the binomial test should be used.")
  dauer<-proc.time()[3] #Measure duration of calculation
  cb<-function(a,b){(a*(C@D[2]+1))+b+1}
  means<-vars<-n<-matrix(C@n,nrow=C@D[1]+1,byrow=TRUE)
  for(i in 1:(C@D[1]+1)){for(j in 1:(C@D[2]+1)){
    means[i,j]<-mean(C@data[[cb(i-1,j-1)]])
    vars[i,j]<-var(C@data[[cb(i-1,j-1)]])
  }}
  padj.fun<-function(x,y1,y2,z){
    rho<-matrix(nrow=C@D[1]*C@D[2],ncol=C@D[1]*C@D[2])
    for(i in 1:C@D[1]){for(j in 1:C@D[2]){for(l in 1:C@D[1]){for(m in 1:C@D[2]){
      if(i==l & j==m) rho[(j-1)*C@D[1]+i,(m-1)*C@D[1]+l]<-1
      if(i==l & j!=m) rho[(j-1)*C@D[1]+i,(m-1)*C@D[1]+l]<-x[i,j]*x[i,m]*sqrt((1-y1[i,j])*(1-y1[i,m]))
      if(i!=l & j==m) rho[(j-1)*C@D[1]+i,(m-1)*C@D[1]+l]<-(1-x[i,j])*(1-x[l,j])*sqrt((1-y2[i,j])*(1-y2[l,j]))
      if(i!=l & j!=m) rho[(j-1)*C@D[1]+i,(m-1)*C@D[1]+l]<-0
    }}}}
    if(C@D[1]==1&&C@D[2]==1) return(1-pmvnorm(lower=rep(-Inf,C@D[1]*C@D[2]),upper=rep(qnorm(1-z),C@D[1]*C@D[2]),mean=rep(0,C@D[1]*C@D[2]),sigma=as.numeric(rho),maxpts=50000,abseps=1E-05))
    else return(1-pmvnorm(lower=rep(-Inf,C@D[1]*C@D[2]),upper=rep(qnorm(1-z),C@D[1]*C@D[2]),mean=rep(0,C@D[1]*C@D[2]),corr=rho,maxpts=50000,abseps=1E-05))
  }
  library(mvtnorm)
  lambda2<-lambda1<-padj<-punadj<-tmin<-matrix(nrow=C@D[1],ncol=C@D[2])
  for(i in 1:C@D[1]){for(j in 1:C@D[2]){
    ta<-(means[i+1,j+1]-means[i+1,1])/sqrt((vars[i+1,j+1]/n[i+1,j+1])+(vars[i+1,1]/n[i+1,1]))
    tb<-(means[i+1,j+1]-means[1,j+1])/sqrt((vars[i+1,j+1]/n[i+1,j+1])+(vars[i+1,1]/n[1,j+1]))
    tmin[i,j]<-min(ta,tb)
    punadj[i,j]<-1-pt(tmin[i,j],df=sum(n)-(C@D[2]+1)*(C@D[1]+1)) #one-sided
    lambda1[i,j]<-n[1,j+1]/(n[i+1,j+1]+n[1,j+1])
    lambda2[i,j]<-n[i+1,1]/(n[i+1,j+1]+n[i+1,1])
  }}
  Omega<-Omega2.fun(C@D[1],C@D[2])
  for(i in 1:C@D[1]){for(j in 1:C@D[2]){
    padj[i,j]<-max(apply(Omega,3,padj.fun,y1=lambda1,y2=lambda2,z=punadj[i,j]))
  }}
  p<-round(as.numeric(t(padj)),4)
  dauer<-proc.time()[3]-dauer
  return(new("mintest",p=p,stat=as.numeric(t(tmin)),gnames=kombi(C@D),method="Hung (2000)",nboot=10000,simerror=.02,duration=dauer))
}

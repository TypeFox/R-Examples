
MaxProLHD <-
function(n,p,s=2,temp0=0,nstarts=1,itermax=400,total_iter=1000000){
  
  m<-n
  k<-p
  
  if(temp0==0){
    avgdist1=1/(n-1)
    avgdist2=(1/((n-1)^(k-1)*(n-2)))^(1/k)
    delta=avgdist2-avgdist1
    temp0=-delta/log(0.99)
  }
  
  t00<-Sys.time()
  aaa<-.C("MaxProLHD",as.integer(m),as.integer(k),as.integer(nstarts),
          as.integer(itermax),as.integer(total_iter),design=integer(m*k),measure=double(1), 
          as.double(temp0),ntotalI=integer(1),as.integer(s),PACKAGE="MaxPro")
  t01<-Sys.time()
  
  time_rec=t01-t00
  
  dd=matrix(aaa$design,ncol=k,nrow=n,byrow=TRUE)
  
  scaled_deisgn<-(apply(dd,2,rank)-0.5)/n
  
  prod_criterion<-function(D)
  {
    logDIST=s*log(dist(D[,1]))
    for(j in 2:p)
      logDIST=logDIST+s*log(dist(D[,j]))
    imin=which.min(logDIST)
    value=-logDIST[imin]+log(sum(exp(logDIST[imin]-logDIST)))
    return(exp((value-log(choose(n,2)))/(p*s)))
  }
  
  measure<-prod_criterion(scaled_deisgn)
  
  val<-list(Design=scaled_deisgn,temp0=temp0,measure=measure,time_rec=time_rec,ntotal=aaa$ntotalI)
  
  return(val)
}

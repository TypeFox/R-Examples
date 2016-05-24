
MaxPro<-function(InitialDesign,s=2,iteration=10){
  
  D0<-as.matrix(InitialDesign)    
  n<-nrow(D0)
  p<-ncol(D0)
  D0<-apply(D0,2, function(x) (x-min(x))/(max(x)-min(x))) 
  
  fgr=function(x)
  {
    D=matrix(x,nrow=n,ncol=p)
    dis=(dist(D[,1]))^s
    for(j in 2:p)
      dis=dis*(dist(D[,j]))^s
    DIST=as.matrix(dis)
    fn=sum(1/dis)
    lfn=log(fn)
    I=diag(n)
    diag(DIST)=rep(1,n)
    A=B=D
    for(j in 1:p)
    {
      A=t(outer(D[,j],D[,j],"-"))
      diag(A)=rep(1,n)
      B[,j]=diag((1/A-I)%*%(1/DIST-I))
    }
    G=s*B/fn
    return(list("objective"=lfn,"gradient"=G))
  }
  
  D1<-D0
  for(i in 1:iteration)
  {
    a=nloptr(c(D1),fgr,opts=list("algorithm"="NLOPT_LD_LBFGS","maxeval"=100),lb=rep(0,n*p),ub=rep(1,n*p))
    D1=matrix(a$sol,nrow=n,ncol=p)
  }
   
  prod_criterion<-function(D)
  {
    logDIST=s*log(dist(D[,1]))
    for(j in 2:p)
      logDIST=logDIST+s*log(dist(D[,j]))
    imin=which.min(logDIST)
    value=-logDIST[imin]+log(sum(exp(logDIST[imin]-logDIST)))
    return(exp((value-log(choose(n,2)))/(p*s)))
  }
  
  Q1=prod_criterion(D1)
  
  val<-list(Design=D1,measure=Q1)
  
  return(val)
  
}




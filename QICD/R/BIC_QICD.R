checkloss<-function(res,tau=0.5) {
  # the check loss function for quantile regression
  Pres=(res+abs(res))/2
  Nres=(res-abs(res))/2
  return(tau*Pres-(1-tau)*Nres)
}


QBIC<-function(y,X,beta,tau=0.5,const=6){
  #QBIC for high dimensional case
  size=dim(X)
  n=size[1]
  p_n=size[2]-1
  #p_n=n
  C_n=log(log(n))/n
  df=sum(abs(beta)>1e-6)
  return(log(sum(checkloss(y-X%*%beta,tau)))+log(p_n)*df*C_n/const)
}

QICD.BIC<-function(y, x, beta=NULL, const=6, tau, lambda, a=3.7,funname="scad",intercept=TRUE,thresh=1e-06,
                  maxin=100,maxout=20,plot.off=F,...)
  #x: input nxp matrix, of dimension nobs x nvars; each row is an observation vector. 
  #y: response variable, length n vector
  #beta: initial value, the defaul value is NULL
  #const is the tuning parameter in QBIC
  #tau is the quantile value
  #lambda is the tuning parameter sequence
  #nfolds: number of folds, default is 10
  #a is scale parameter, the default value is 3.7 for SCAD
  #funname is the name of nonconvex penalty function, could be scad, mcp and lasso, the default
  #value is scad
  #intercept is a logical value,should intercept(s) be fitted (default=TRUE) or set to zero(FALSE)
  #thresh is the convergence threshold for coordinate descent and majorization minimization step.
#Default value is 1E-6 
#maxin: maximum number of iterations for inside coordinate descent,default value is 100
#maxout: maximum number of iterations for outside MM step,default value is 20
#plot.off: control a plot.Default is False, if TRUE, a plot will be given
#...: arguments could be passed to plot
{
  HBIC=NULL
  nzero=NULL
  nlambda=length(lambda)
  n=dim(x)[1]
  for (j in 1:nlambda){
    res=QICD(y,x,tau=tau,beta=beta,lambda=lambda[j],a=a,funname=funname,
                   intercept=intercept,thresh=thresh,maxin=maxin,maxout=maxout)
    beta=res$beta_final
    if (intercept)
      HBIC=c(HBIC,QBIC(y,cbind(x,rep(1,n)),beta,tau,const))
    else
      HBIC=c(HBIC,QBIC(y,x,beta,tau,const))
    nzero=c(nzero,res$df)
  }
  HBIC.min=min(HBIC)
  #HBICsd=sd(HBIC)
  #HBICup=HBIC+HBICsd
  #HBIClo=HBIC-HBICsd  
  lambda.min=lambda[which.min(HBIC)]
  #lambda.1se=max(lambda[which(abs(HBIC-HBIC.min)<HBICsd)])
  if (!plot.off){
    plot(lambda,HBIC,type="l",...)
    #lines(lambda,HBICup,lty=3)
    #lines(lambda,HBIClo,lty=3)
  }  
  obj=list (lambda=lambda,HBIC=HBIC,nzero=nzero,
        lambda.min=lambda.min)
  class(obj)="QICD.BIC"
  obj
}
#source("QICD.R")
checkloss<-function(res,tau=0.5) {
  # the check loss function for quantile regression
  Pres=(res+abs(res))/2
  Nres=(res-abs(res))/2
  return(tau*Pres-(1-tau)*Nres)
}

QICD.cv<-function(y, x, beta=NULL, tau, lambda,nfolds=10, a=3.7,funname="scad",intercept=TRUE,thresh=1e-06,
               maxin=100,maxout=20,mc.cores=getOption("mc.cores", 1L),plot.off=F,...)
#x: input nxp matrix, of dimension nobs x nvars; each row is an observation vector. 
#y: response variable, length n vector
#beta: initial value, the defaul value is NULL
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
#mc.cores: The number of cores to use for parallel computing.
#plot.off: control a plot.Default is False, if TRUE, a plot will be given
#...: arguments could be passed to plot
{
  size=dim(x)
  n=size[1]
  p=size[2]
  nlambda=length(lambda)
  cross_index=sample(1:nfolds,size=n,replace=TRUE)
  index=sample(1:n,n)
  cvm=NULL
  cvsd=NULL
  nzero=NULL
  for(j in 1:nlambda){
    prediction_test=NULL
    nzero_test=NULL
    
    cv_QICD<-function(fold){
        x_train=x[cross_index!=fold,]
        y_train=y[cross_index!=fold]
        x_test=x[cross_index==fold,]
        y_test=y[cross_index==fold]
        res_train=QICD(y_train,x_train,tau=tau,beta=beta,lambda=lambda[j],a=a,funname=funname,
                        intercept=intercept,thresh=thresh,maxin=maxin,maxout=maxout)
        return(res_train)    
    }
    cv_obj=mclapply(1:nfolds,cv_QICD,mc.cores = mc.cores)
    for (i in 1:nfolds){
        x_test=x[cross_index==i,]
        y_test=y[cross_index==i]
        res_train=cv_obj[[i]]
        beta_train=res_train$beta_final
        nbeta=dim(x_test)[1]
        #training coeffcients 
        if (intercept)
          prediction_test=c(prediction_test,sum(checkloss(y_test-cbind(x_test,rep(1,nbeta))%*%beta_train,tau)))
        else
          prediction_test=c(prediction_test,sum(checkloss(y_test-x_test%*%beta_train,tau)))
        nzero_test=c(nzero_test,res_train$df)
    }
    ave_prediction=mean(prediction_test)
    ave_nzero=mean(nzero_test)
    cvm=c(cvm,ave_prediction)
    cvsd=c(cvsd,sd(prediction_test)/sqrt(nfolds))
    nzero=c(nzero,ave_nzero)
  }
  cvm.min=min(cvm)
  #cvsd=sd(cvm)
  cvup=cvm+cvsd
  cvlo=cvm-cvsd
  lambda.min=lambda[which.min(cvm)]
  lambda.1se=max(lambda[which(abs(cvm-cvm.min)<cvsd[which.min(cvm)])])
  if (!plot.off){
    plot.args=list(x=lambda,
                   y=cvm,ylim=range(cvup,cvlo),
                   xlab="lambda",ylab="cvm",type="l")
    lines.cvup.args=list(x=lambda,
                         y=cvup,lty=3)
    lines.cvlo.args=list(x=lambda,
                         y=cvlo,lty=3)
    # plot(lambda,cvm,type="l",...)
    # lines(lambda,cvup,lty=3,...)
    # lines(lambda,cvlo,lty=3,...)
    new.args=list(...)
    if(length(new.args)){
      plot.args[names(new.args)]=new.args
      lines.cvup.args[names(new.args)]=new.args
      lines.cvlo.args[names(new.args)]=new.args
    }
    do.call("plot",plot.args)
    do.call("lines",lines.cvup.args)
    do.call("lines",lines.cvlo.args)
  }
  
  obj=list (lambda=lambda,cvm=cvm,cvsd=cvsd,cvup=cvup,cvlo=cvlo,nzero=nzero,lambda.min=lambda.min,lambda.1se=lambda.1se)
  class(obj)="QICD.cv"
  obj
}
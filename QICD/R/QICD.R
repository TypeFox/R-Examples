#.First.lib <- function(lib, pkg)
#{
#  library.dynam("LIBNAME", pkg, lib)
#}
QICD<-function(y, x, beta=NULL, tau, lambda, a=3.7,funname="scad",intercept=TRUE,thresh=1e-06,
               exclude=NULL,maxin=100,maxout=20)
#x: input nxp matrix, of dimension nobs x nvars; each row is an observation vector. 
#y: response variable, length n vector
#beta: initial value, the defaul value is NULL
#tau is the quantile value
#lambda is the tuning parameter sequence
#a is scale parameter, the default value is 3.7 for SCAD, no a is applicable if penalty is lasso
#funname is the name of nonconvex penalty function, could be scad, mcp and lasso, the default
#value is scad
#intercept is a logical value,should intercept(s) be fitted (default=TRUE) or set to zero(FALSE)
#thresh is the convergence threshold for coordinate descent and majorization minimization step.
#Default value is 1E-6
#exclude: indices of variables to be excluded from the model  
#maxin: maximum number of iterations for inside coordinate descent,default value is 100
#maxout: maximum number of iterations for outside MM step,default value is 20
{
  #dyn.load("QCD.dll")
  p=ncol(x)
  #raw dimension of data matrix x
  exclude_zero=apply(x,2,allzero)
  if (is.logical(exclude)){
    exclude=(exclude|exclude_zero)
    x=x[,!exclude]
  }
  else{
    exclude_temp=rep(F,p)
    exclude_temp[exclude]=T
    exclude=(exclude_temp|exclude_zero)
    x=x[,!exclude]
  } 
    
#######################################################
  if (length(x)==0)
    stop("x is an all zeros matrix")
  
  nyrow<-as.integer(length(y))
  nxcol<-as.integer(ncol(x))
  if(funname=="scad")
	  index<-0  
  else if (funname=="mcp")
	  index<-1
  else if (funname=="lasso")
    index<-2
  else
    stop("wrong penalty function")
  if (intercept){
    p=p+1
    x=cbind(x,rep(1,nyrow))
    #create observation matrix with the last column to be ones
    index1=1
    #intercept indicator
    if (is.null(beta))
      beta=rep(0,nxcol+1)
    else
      beta=beta[!exclude]
  }
  else{
    index1=0
    if (is.null(beta))
      beta=rep(0,nxcol)
    else
      beta=beta[!exclude]
  }
  index<-as.integer(index)
  lambda=lambda/nyrow    
  nlambda=length(lambda)
  beta_final=NULL
  df=NULL
  #none zero numbers for coefficients
  for (j in 1:nlambda){
    beta1=beta
    i=0
    repeat{
      beta0<-beta1
      pre_value<-x%*%beta1
      out<-.C("QCD",as.double(y),as.double(x),as.double(beta0),as.double(beta1)
              ,as.double(pre_value),nyrow,nxcol,as.double(tau),as.double(lambda[j]),
              as.double(a),index,as.integer(index1),as.double(thresh),as.integer(maxin))
      beta1<-out[[4]]
      i<-i+1
      distance=sqrt(sum((beta1-beta0)^2))
      if ((i>maxout)| (distance<thresh)) break
    }
    df=c(df,sum(abs(beta1)>thresh))
    beta_temp=rep(0,p)
    beta_temp[!exclude]=beta1
    beta_final=cbind(beta_final,beta_temp)
    #final beta coefficients
  }
  obj=list (beta_final=beta_final,lambda=lambda,df=df,dim=dim(beta_final))
  class(obj)="QICD"
  obj
}

allzero<-function(x)
# are all x column values are zeros
{
  return(all(x==0))
}



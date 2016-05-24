#----------------------------------------------------------------------------------#
# Package: camel                                                                   #
# camel.slim(): The user interface for slim()                                      #
# Author: Xingguo Li                                                               #
# Email: <xingguo.leo@gmail.com>                                                   #
# Date: Aug 23th, 2013                                                             #
# Version: 0.1.0                                                                   #
#----------------------------------------------------------------------------------#

camel.slim <- function(X, 
                       Y, 
                       lambda = NULL,
                       nlambda = NULL,
                       lambda.min.ratio = NULL,
                       method="lq",
                       q = 2,
                       prec = 1e-4,
                       max.ite = 1e4,
                       mu = 0.01,
                       intercept = TRUE,
                       verbose = TRUE)
{
  if(method!="dantzig" && method!="lq"){
    cat("\"method\" must be dantzig or lq lasso.\n")
    return(NULL)
  }
  if(method=="lq"){
    if(q!=1 && q!=2){
      cat("q must be either 1 or 2 when method = \"lq\".\n")
      return(NULL)
    }
  } else q=0
  if(verbose) {
    cat("Sparse Linear Regression with L1 Regularization.\n")
  }
  n = nrow(X)
  d = ncol(X)
  maxdf = max(n,d)
  
  xm=matrix(rep(.colMeans(X,n,d),n),nrow=n,ncol=d,byrow=T)
  x1=X-xm
  sdx=sqrt(diag(t(x1)%*%x1)/(n-1))
  Cxinv=diag(1/sdx)
  xx=x1%*%Cxinv
  ym=mean(Y)
  y1=Y-ym
  sdy=sqrt(sum(y1^2)/(n-1))
  yy=y1/sdy
  corr=cor(xx,yy)
  
  if(intercept){
    xx = cbind(rep(1, nrow(xx)), xx)
    d = d+1
  }
  
  if(!is.null(lambda)) nlambda = length(lambda)
  if(is.null(lambda)){
    if(is.null(nlambda))
      nlambda = 5
    if(is.null(lambda.min.ratio)){
      if(method=="dantzig")
        lambda.min.ratio = 0.8
      else
        lambda.min.ratio = 0.25
    }
    if(method=="dantzig")
      lambda.max = max(corr)
    else
      lambda.max = pi*sqrt(log(d)/n)
      
    lambda.min = lambda.min.ratio*lambda.max
    lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
    rm(lambda.max,lambda.min,lambda.min.ratio)
    gc()
  }
  begt=Sys.time()
  if(method=="dantzig") # dantzig
    out = camel.slim.dantzig.mfista(yy, xx, lambda, nlambda, n, d, maxdf, mu, max.ite, prec,intercept,verbose)

  if(method=="lq") {
    if(is.null(q)) q=2;
    if(q==1) { # lad lasso
      out = camel.slim.lad.mfista(Y, xx, lambda*n, nlambda, n, d, maxdf, mu, max.ite, prec,intercept,verbose)
    }
    if(q==2) { # sqrt lasso
      out = camel.slim.sqrt.mfista(Y, xx, lambda*sqrt(n), nlambda, n, d, maxdf, mu, max.ite, prec,intercept,verbose)
    }
  }
  runt=Sys.time()-begt
  
  sparsity=rep(0,nlambda)
  for(i in 1:nlambda)
    sparsity[i] = sum(out$beta[[i]]!=0)/d
  
  est = list()    
  intcpt=matrix(0,nrow=1,ncol=nlambda)
  
  if(!intercept){
    beta1=matrix(0,nrow=d,ncol=nlambda)
    if(q==1 || q==2){
      for(k in 1:nlambda){
        tmp.beta = out$beta[[k]]
        intcpt[k]=-xm[1,]%*%Cxinv%*%tmp.beta
        beta1[,k]=Cxinv%*%tmp.beta
      }
    }
    else{
      for(k in 1:nlambda){
        intcpt[k]=ym-xm[1,]%*%Cxinv%*%out$beta[[k]]*sdy
        beta1[,k]=Cxinv%*%out$beta[[k]]*sdy
      }
    }
  } else {
    beta1=matrix(0,nrow=d-1,ncol=nlambda)
    if(q==1 || q==2){
      for(k in 1:nlambda){
        tmp.beta = out$beta[[k]][2:d]
        intcpt[k]=out$beta[[k]][1]-xm[1,]%*%Cxinv%*%tmp.beta
        beta1[,k]=Cxinv%*%tmp.beta
      }
    }
    else{
      for(k in 1:nlambda){
        tmp.beta = out$beta[[k]][2:d]
        intcpt[k]=ym-xm[1,]%*%Cxinv%*%tmp.beta*sdy+out$beta[[k]][1]*sdy
        beta1[,k]=Cxinv%*%tmp.beta*sdy
      }
    }
  }
  
  est$beta = beta1
  est$intercept = intcpt
  est$Y = Y
  est$X = X
  est$lambda = lambda
  est$nlambda = nlambda
  est$sparsity = sparsity
  est$method = method
  est$q = q
  est$ite =out$ite
  est$verbose = verbose
  est$runtime = runt
  class(est) = "slim"
  return(est)
}

print.slim <- function(x, ...)
{  
  cat("\n slim options summary: \n")
  cat(x$nlambda, " lambdas used:\n")
  print(signif(x$lambda,digits=3))
  cat("Method=", x$method, "\n")
  if(x$method=="lq"){
    if(x$q==1){
      cat("q=",x$q," loss, LAD Lasso\n")
    } else {
      if(x$q==2)
        cat("q=",x$q," loss, SQRT Lasso\n")
      else
        cat("q=",x$q," loss\n")
    }
  }
  cat("Sparsity level:",min(x$sparsity),"----->",max(x$sparsity),"\n")
  if(units.difftime(x$runtime)=="secs") unit="secs"
  if(units.difftime(x$runtime)=="mins") unit="mins"
  if(units.difftime(x$runtime)=="hours") unit="hours"
  cat("Runtime:",x$runtime," ",unit,"\n")
}

plot.slim <- function(x, ...)
{
  matplot(x$lambda, t(x$beta), type="l", main="Regularization Path",
          xlab="Regularization Parameter", ylab="Coefficient")
}

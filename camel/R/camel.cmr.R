#----------------------------------------------------------------------------------#
# Package: camel                                                                   #
# camel.cmr(): The user interface for crm()                                        #
# Author: Xingguo Li                                                               #
# Email: <xingguo.leo@gmail.com>                                                   #
# Date: Seq 3rd, 2013                                                              #
# Version: 0.2.0                                                                   #
#----------------------------------------------------------------------------------#

camel.cmr <- function(X, 
                      Y, 
                      lambda = NULL,
                      nlambda = NULL,
                      prec = 1e-3,
                      max.ite = 1e3,
                      mu = 0.01,
                      verbose = TRUE)
{
  if(verbose) {
    cat("Multivariate Regression with Calibration via MFISTA.\n")
  }
  n = nrow(X)
  d = ncol(X)
  m = ncol(Y)
  maxdf = max(n,d)
  
  xm=matrix(rep(.colMeans(X,n,d),n),nrow=n,ncol=d,byrow=T)
  x1=X-xm
  sdx=sqrt(diag(t(x1)%*%x1)/(n-1))
  Cxinv=diag(1/sdx)
  xx=x1%*%Cxinv
  ym=matrix(rep(.colMeans(Y,n,m),n),nrow=n,ncol=m,byrow=T)
  yy=Y-ym
#   xx = X
#   yy = Y
  
  if(!is.null(lambda)) nlambda = length(lambda)
  if(is.null(lambda)){
    if(is.null(nlambda))
      nlambda = 10
    lambda.max = sqrt(log(d))+sqrt(m)
    
    lambda = 2^(seq(5, -4, length = nlambda))*lambda.max
    rm(lambda.max)
    gc()
  }
  begt=Sys.time()
  out = camel.cmr.mfista(yy, xx, lambda, nlambda, n, d, m, mu, max.ite, prec)

  runt=Sys.time()-begt
  
  sparsity=rep(0,nlambda)
  for(i in 1:nlambda)
    sparsity[i] = sum(out$beta[[i]]!=0)/(d*m)
  
  est = list()    
  
  beta1 = vector("list", nlambda)
  intcpt = vector("list", nlambda)
  for(k in 1:nlambda){
    tmp.beta = out$beta[[k]]
    intcpt[[k]]=ym[1,]-xm[1,]%*%Cxinv%*%tmp.beta
    beta1[[k]]=Cxinv%*%tmp.beta
  }
  
#   est$beta = out$beta
  est$beta = beta1
  est$Y = Y
  est$X = X
  est$lambda = lambda
  est$nlambda = nlambda
  est$sparsity = sparsity
#   est$method = method
#   est$q = q
  est$ite =out$ite
  est$verbose = verbose
  est$runtime = runt
  class(est) = "cmr"
  return(est)
}

print.cmr <- function(x, ...)
{  
  cat("\n cmr options summary: \n")
  cat(x$nlambda, " lambdas used:\n")
  print(signif(x$lambda,digits=3))
  cat("Sparsity level:",min(x$sparsity),"----->",max(x$sparsity),"\n")
  if(units.difftime(x$runtime)=="secs") unit="secs"
  if(units.difftime(x$runtime)=="mins") unit="mins"
  if(units.difftime(x$runtime)=="hours") unit="hours"
  cat("Runtime:",x$runtime," ",unit,"\n")
}


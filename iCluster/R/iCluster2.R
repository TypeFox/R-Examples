###############################################################
#iCluster   version 2.1.0 lasso with variance weighted shrinkage
#last update: 05/01/2012
#
#ARGUMENTS:
#datasets: multiple genomic platform (MGP) data in list format.
#Arrange each data set such that rows are samples and columns are genes.
#k: number of clusters.
#lambda: a list of m vectors containing penalty parameters
#max.iter: maximum number of iteration in the EM algorithm.
#scalar: a logical argument specifying whether scalar covariance matrix should be used
#VALUE:
#expZ: MLE estimate of E[Z*|X] of dimension K-1 by n where n is sample size.
#clusters: cluster assignments.
################################################################

iCluster2 = function(datasets, k, lambda=NULL, scale=T, scalar=F, max.iter=10, verbose=T){

  n=nrow(datasets[[1]])
  m=length(datasets)
  p=unlist(lapply(1:m,function(l){ncol(datasets[[l]])}))
  sum.p=sum(p)
    
  #pre-checks
  if(is.list(datasets)==F)stop("Datasets must be a list")
  if(sum(is.na(datasets))>0)stop("Data cannot have NAs. please exclude or impute missing values.")
  if(missing(k))stop("Must specify the number of clusters k")

  #create a big matrix combining multiple data types
  stacked=do.call(cbind, datasets)
  if(scale){stacked=scale(stacked,center=T,scale=F)}
  
  if(verbose)cat("Generating warm start...",'\n')
  
  #initialize W_sum.px(k-1) and PSI_sum.pxsum.p using PCA solution
  svd.x=svd(stacked) 
  W=svd.x$v[,1:(k-1),drop=F]%*%diag(sqrt(svd.x$d[1:(k-1)]),nrow=(k-1),ncol=(k-1)) 
  d=svd.x$d^2
  sigma2=sum(svd.x$d[-c(1:k-1)])/(n-k+1)   

  PSI=rep(sigma2,sum.p)
  #PSI=rep(1,sum.p)

  sigma=W%*%t(W)
  diag.idx=1+0:(sum.p-1)*(sum.p+1)
  sigma[diag.idx]=sigma[diag.idx]+PSI
  inv.sigma=ginv(sigma)
  #inv.sigma=solve(sigma)
  

 #setting default values if penalty type and/or lambda not specified
  if(is.null(lambda)){lambda=as.list(rep(0.5,m))}

  expZ=matrix(rnorm((k-1)*n),nrow=k-1,ncol=n)
  clusters=rep(1,n)
  
  iter=1
  conv.rate=1
  save.conv.rate=NULL
  RI=NULL
  RI.iter=0
  if(verbose)cat(paste("K=",k,":",sep=""))
  while((RI.iter<0.95)&(iter<=max.iter)){
  #while(iter<max.iter){
  if(verbose)cat(iter)

  #E-step
  expZ.old=expZ
  expZ=t(W)%*%inv.sigma%*%t(stacked)   #kxn
  expZ=multi.l2n(expZ,row.norm=T)
  varZ=diag(c(1),nrow=k-1,ncol=k-1)-t(W)%*%inv.sigma%*%W
  expZZ=varZ+expZ%*%t(expZ)
  inv.expZZ=solve(expZZ)

  #M-step
  W=t(stacked)%*%t(expZ)%*%inv.expZZ
  for(i in 1:m){
      if(i==1){idx=1:p[i]}else{idx=(sum(p[1:(i-1)])+1):sum(p[1:i])}
      aa=PSI[idx]
          	for(j in 1:(k-1)){
            	lam=BinarySearch(W[idx,j,drop=F]*2/aa, lambda[[i]]*sqrt(p[i]))
            	lam=lam*aa/2
            	W[idx,j]=soft(W[idx,j,drop=F], lam)
          	}
	  W[idx,]=multi.l2n(W[idx,,drop=F],row.norm=F)
	
  }
  resid=t(stacked)-W%*%expZ
  a=diag(resid%*%t(resid))
  if(scalar){PSI=rep(sum(a)/n/sum.p,length(a))}else{
  PSI=a/n
  }

  sigma=W%*%t(W)
  sigma[diag.idx]=sigma[diag.idx]+PSI
  inv.sigma=ginv(sigma)
  #inv.sigma=solve(sigma)

  #monitor convergence based on expZ
  conv.rate=max(abs(expZ-expZ.old))
  save.conv.rate[iter]=conv.rate
  
  #monitor convergence based on cluster ass
  clusters.old=clusters
  clusters=kmeans(t(expZ),k,nstart=100)$cluster
  RI.iter=RandIndex(clusters.old, clusters)
  RI[iter]=RI.iter

  iter=iter+1
  }
  if(verbose)cat('\n')
  if(iter==max.iter){warning("Algorithm didn't converge. Check convergence history fit$RI. Cluster assignments may not be stable. Try increase the number of EM iterations by max.iter")}
  kmeans.fit=kmeans(t(expZ),k,nstart=100)
  output=list(expZ=expZ, W=W, PSI=PSI, clusters=kmeans.fit$cluster,centers=kmeans.fit$centers, RI=RI, lambda=lambda)
  return(output)
}

###soft thresholding###
soft=function (a, para){
  b <- sort(abs(a))
  b <- abs(a) - para
  b <- (b + abs(b))/2 #|b|+ take the positive value
  b <- sign(a) * b
  b
}

multi.l2n=function(a, row.norm=T){
if(row.norm){norma=sqrt(apply(a,1,function(x){sum(x^2)}))}else{
norma=sqrt(apply(a,2,function(x){sum(x^2)}))}
norma[norma==0]=.05
a=t(t(a)/norma)
return(a)
}

l2n <- function(vec){
  a <- sqrt(sum(vec^2))
  if(a==0) a <- .05
  return(a)
}

RandIndex=function (c1, c2)
{
    c1 <- as.vector(c1)
    c2 <- as.vector(c2)
    xx <- outer(c1, c1, "==")
    yy <- outer(c2, c2, "==")
    upper <- row(xx) < col(xx)
    xx <- xx[upper]
    yy <- yy[upper]
    a <- sum(as.numeric(xx & yy))
    d <- sum(as.numeric(!xx & !yy))
    (a+d)/choose(length(c2),2)
}
BinarySearch <- function(argu,sumabs){
  if(l2n(argu)==0 || sum(abs(argu/l2n(argu)))<=sumabs) return(0)
  lam1 <- 0
  lam2 <- max(abs(argu))-1e-5
  iter <- 1
  while(iter < 150){
    su <- soft(argu,(lam1+lam2)/2)
    if(sum(abs(su/l2n(su)))<sumabs){
      lam2 <- (lam1+lam2)/2
    } else {
      lam1 <- (lam1+lam2)/2
    }
    if((lam2-lam1)<1e-6) return((lam1+lam2)/2)
    iter <- iter+1
  }
  warning("Didn't quite converge")
  return((lam1+lam2)/2)
}

msqrt <- function(x){
  eigenx <- eigen(x)
  return(eigenx$vectors%*%diag(sqrt(pmax(0,eigenx$values)))%*%t(eigenx$vectors))
}

fastsvd <- function(x,z){
  # fast svd of t(x)%*%z, where ncol(x)>>nrow(x) and same for z
  xx=x%*%t(x)
  xx2=msqrt(xx)
  y=t(z)%*%xx2
  a=svd(y)
  v=a$u
  d=a$d
  zz=z%*%t(z)
  zz2=msqrt(zz)
  y=t(x)%*%zz2
  a=svd(y)
  u=a$u
  return(list(u=u,v=v,d=d))
}


ginv=function(X, tol = sqrt(.Machine$double.eps))
{
  ## Generalized Inverse of a Matrix
  dnx <- dimnames(X)
  if(is.null(dnx)) dnx <- vector("list", 2)
  s <- svd(X)
  nz <- s$d > tol * s$d[1]
  structure(
    if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else X,
    dimnames = dnx[2:1])
}

RHT.fun <-
function(path.idx, dat, nsim=1000, seed=123){


my.HotellingT <-
function(X, u=0, lambda=1)
{


  my.InvCov<-function(X, lambda=1)
  {
    p=ncol(X)
    N=nrow(X)
    n=N-1
    ##### covariance
    X.cov=cov(X, X, use="pairwise.complete.obs")
    diag(X.cov)[is.na(diag(X.cov))] <- 1
    X.cov[is.na(X.cov)]=0
  
    svd.x = svd(X.cov)
    r = sum(svd.x$d>10^(-6))
    Dinv <- matrix(0,r,r)
    diag(Dinv)=1/(svd.x$d[1:r] +lambda)
    H <- as.matrix(svd.x$u[, 1:r])
    X.InvCov <- H%*%Dinv%*%t(H)
    return(X.InvCov)
  }



  temp=apply(!is.na(X), 2, sum)
  X.m=X[,temp>=1]

  p=ncol(X.m)
  N=nrow(X.m)
  n=N-1
  col.mean=colMeans(X.m, na.rm=TRUE)

  ##### sd
  W=my.InvCov(X.m, lambda=lambda)

  result=N*t(col.mean-u)%*%W%*%(col.mean-u)
  return(result)
  
}



get.lambda <-
function(d, use.norm=NULL){

  get.moment <- function(d, lambda){
    const <- d/(d+lambda)
    p <- length(d)/2
    m1 <- sum(const)
    m2 <- sum(const^2)*2+m1^2
    return(c(m1,m2))
  }
 
  ff <- function(r, v, R2){
    abs(R2-gamma(v)*gamma(2*r+v)/gamma(r+v)^2)
  }

  get.para <- function(d, lambda){
    m <- get.moment(d, lambda)
    R2=m[2]/m[1]^2
    p <- length(d)
    r <- optimize(ff, interval=c(0, 5), v=p/2, R2=R2, maximum=FALSE)$minimum
    A <- m[1]*gamma(p/2)/gamma(r+p/2)/2^r
    return(c(A, r, p))
  }

  gg <- function(lambda, d, alpha=0.95){
    para <- get.para(d, lambda)
    A <- para[1]; r <- para[2]; p <- para[3]
    1/(1+lambda)*qchisq(alpha, p) - A*qchisq(alpha,p)^r
  }

  ggl <- function(lambda, d, alpha=0.95){
    p=length(d)
    qnull= qnorm(alpha, p/(1+lambda), sqrt(2*p)/(1+lambda))
    obs.m= sum(d/(d+lambda))
    obs.v = sum(2* (d/(d+lambda))^2)
    qobs=qnorm(alpha, obs.m, sqrt(obs.v))
    return(qnull-qobs)
  }

  if (is.list(d)){
      d=lapply(d, function(di){
        di <- di/sum(di)*length(di)
      })
      d=unlist(d)
  } else {
       d <- d/sum(d)*length(d)
  }
  val <- NULL
  lamb.r <- (1:10000)/100
  if (length(d)>=100) use.norm=TRUE else use.norm=FALSE
  for (lambda in lamb.r){
       if (use.norm){
         val <- c(val, ggl(lambda,d))   
       } else {
         val <- c(val, gg(lambda,d))
       }  
  }
  i1 <- which.max(val)
  lambda.c <- lamb.r[i1]
  lambda <- lambda.c*sum(d)/length(d)
  return(lambda)
}




    set.seed(seed)
    nP <- length(path.idx)
    pval = rep(NA, nP)
    lambs = rep(0, nP)
    
    
    for (i in 1:nP){
      path.i <- path.idx[[i]]
      
      X <-  dat[, path.i]
    
      nnas <- apply(X, 2, function(x) sum(!is.na(x)))
    
      if (sum(nnas>=2)>=2){
        X=as.matrix(X[, nnas>=1])
        p = ncol(X)     
        n = nrow(X)-1       
	if (p<n) {
	   boot="no" 
	   Nset=1
	   splitSet=list()
           splitSet[[1]]=1:p
	} else {
           Nset=as.integer(5*p/(n-1))
           splitSet = lapply(1:Nset, function(x) sample(1:p, n-1, replace=FALSE))
        }

        d.list = lapply(1:Nset, function(j){
              idxi = splitSet[[j]]
              pi = length(idxi)
              Xi = X[,idxi]

              Xi.cov=cov(Xi, Xi, use="pairwise.complete.obs")
              diag(Xi.cov)[is.na(diag(Xi.cov))] <- 1
              Xi.cov[is.na(Xi.cov)]=0
                 
              svd.x = svd(Xi.cov)
              r=sum(svd.x$d>=10^(-6))
              return(svd.x$d[1:r]) })
        
        
        if (Nset>1) use.norm=TRUE else use.norm=FALSE
        lamb=get.lambda(d.list, use.norm=use.norm)
        lambs[i] <- lamb
        
        stat <- 0
        for (j in 1:Nset){
           idxi = splitSet[[j]]
           Xi = X[,idxi]
           stat =stat +my.HotellingT(Xi, 0, lambda=lamb)        
        } 
        
    
        stat0 = rep(0, nsim)
        for (s in 1:nsim){
          X0 <- X
          X0[!is.na(X0)] <- sample(dat[!is.na(dat)], sum(!is.na(X0)))
         
          for (j in 1:Nset){
                 idxi = splitSet[[j]]
                 pi = length(idxi)
                 Xi = X0[,idxi]
                 stat0[s]=stat0[s]+my.HotellingT(Xi, 0, lambda=lamb)         
          } 
        }
             
    pval[i] <- mean(stat0>=stat[1])          
  }
  }
  names(pval) <- names(path.idx)
  
  return(pval)
}


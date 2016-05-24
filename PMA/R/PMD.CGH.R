ChooseLambda1Lambda2 <- function(y,chrom=NULL){ # using as of 8/25/09
  if(is.null(chrom)) chrom <- rep(1, length(y))
  y.lo <- lowess(y, f=min(50/length(y),.2))$y
  y.lo <- soft(y.lo, quantile(abs(y.lo), .45))
  s1=sum(abs(y.lo)) 
  lambdas=exp(seq(log(10^(-4)), .2, len=50))
  out <- matrix(NA, nrow=length(y),ncol=length(lambdas))
  for(chr in unique(chrom)){
    out[chrom==chr,] <- t(FLSA(y[chrom==chr],lambda1=0,lambda2=lambdas))
  }
  for(i in 1:length(lambdas)) out[,i] <- soft(out[,i], lambdas[i])
  s1s <- apply(abs(out),2,sum)
  errs <- abs(s1s-s1)
  bestlam <- lambdas[which.min(errs)]
  return(bestlam)
}




ChooseLambda1Lambda2old <- function(y,chrom=NULL){ # old as of 8/25/09
  if(is.null(chrom)) chrom <- rep(1, length(y))
  y.lo <- lowess(y, f=min(50/length(y),.2))$y
  y.lo <- soft(y.lo, quantile(abs(y.lo), .45))
  s1=sum(abs(y.lo)) 
  lambdas=exp(seq(log(10^(-4)), .2, len=20))
  out <- matrix(NA, nrow=length(y),ncol=length(lambdas))
  for(i in 1:length(lambdas)){
    for(chr in unique(chrom)){
      out[chrom==chr,i] <- FLSA(y[chrom==chr],lambda1=lambdas[i],lambda2=lambdas[i])[1,1,]
#      flsa.out <- diag.fused.lasso.new(y[chrom==chr], lam1=lambdas[i])
#      lam2ind <- which.min(abs(lambdas[i]-flsa.out$lam2))
#      out[chrom==chr,] <- flsa.out$coef[,lam2ind]
    }
  }
  s1s <- apply(abs(out),2,sum)
  errs <- abs(s1s-s1)
  bestlam <- lambdas[which.min(errs)]
  return(bestlam)
}


CGH.SMD <- function(x,lam1,lam2,sumabsu,chrom,niter=20,v, trace, upos, uneg){
  nas <- is.na(x)
  xoo <- x
  xoo[nas] <- mean.na(x)
  vold <- rnorm(length(v))
  for(i in 1:niter){
    if(sum(abs(vold-v))>1e-5 && sum(abs(v))!=0){
      if(trace) cat(i,fill=F)
      argu <- xoo%*%v
      if(upos) argu <- pmax(argu,0)
      if(uneg) argu <- pmin(argu,0)
      lamu <- BinarySearch(argu,sumabsu)
      u <- soft(argu,lamu)/l2n(soft(argu,lamu))
      vnew <- numeric(ncol(x))
      for(j in sort(unique(chrom))){
        xoou <- as.numeric(t(xoo[,chrom==j])%*%u)
        coefs <- FLSA(xoou/l2n(xoou),lam1,lam2)
#        flsa.out <- diag.fused.lasso.new(xoou/l2n(xoou), lam1=lam1)
#        lam2ind <- which.min(abs(flsa.out$lam2-lam2))
#        coefs <- flsa.out$coef[,lam2ind]
        vnew[chrom==j] <- coefs
      }
      vold <- v
      v <- vnew/l2n(vnew) # Doing this re-scaling s.t. v retains norm 1
      if(sum(is.na(v))>0) v <- rep(0, ncol(x))
    }
  }
  if(sum(abs(v))==0) u <- matrix(0, nrow=nrow(x),ncol=1)
  d <- as.numeric(matrix(u,nrow=1)%*%xoo%*%matrix(v, ncol=1))
  return(list(u=u,v=v,d=d))
}

MultiCGH.SMD <- function(x,K,lam1,lam2,sumabsu,chrom,niter,v, trace, upos, uneg){
  uans = vans = dans = NULL
  xres <- x
  for(k in 1:K){
    if(K==1) out <- CGH.SMD(xres,lam1,lam2,sumabsu,chrom,niter,v, trace=trace, upos=upos, uneg=uneg)
    if(K>1) out <- CGH.SMD(xres,lam1,lam2,sumabsu,chrom,niter,v[,k], trace=trace, upos=upos, uneg=uneg)
    uans <- cbind(uans,out$u)
    vans <- cbind(vans,out$v)
    dans <- c(dans, out$d)
    xres <- xres - out$u%*%out$d%*%t(out$v)
  }
  return(list(u=uans,v=vans, d=dans, v.init=v, lambda=lam1, sumabsu=sumabsu, K=K))
}

PMDL1FL <- function(x,K=1,sumabsu=5,lambda=NULL,chrom=NULL,niter=20, v=NULL, trace=TRUE, center=TRUE, rnames=NULL, cnames=NULL, upos, uneg){
  call <- match.call()
  if(center){
    meanx <- mean.na(x)
    x <- x-meanx
  } else {
    meanx <- NULL
  }
  if(is.null(chrom)) chrom <- rep(1, ncol(x))
  if(sumabsu<1 || sumabsu>sqrt(nrow(x))) stop("sumabsu must be at least 1, and no greater than sqrt(nrow(x))")
  if(!is.null(lambda) && lambda<0) stop("lambda must be non-negative")
  v <- CheckPMDV(v,x,K)
  if(is.null(lambda)) lambda <- ChooseLambda1Lambda2(as.numeric(v[,1])) 
  out <- MultiCGH.SMD(x,K=K,lam1=lambda,lam2=lambda,sumabsu=sumabsu,chrom=chrom,niter=niter, v=v, trace=trace, upos=upos, uneg=uneg)
  out$rnames <- rnames
  out$cnames <- cnames
  out$call <- call
  out$meanx <- meanx
  class(out) <- "PMDL1FL"
  if(trace) cat(fill=TRUE)
  return(out)
}

#Note this R code fits the hierarchical normal-normal model using ADM
#Author: Joseph Kelly

gr<-function(y,se,X,mu,confidence.lvl=0.95,intercept=T,eps=0.0001, normal.CI = FALSE){

  ##define some values that will be used often
  muknown <- !missing(mu)
  if (!muknown) {
    prior.mean <- NA
  } else {
	prior.mean <- mu
  }
  type <- 1
  k <- length(y)
  V <- se^2
  if(missing(X)){
    X.ini<-NA
    X <- as.matrix(rep(1,k))
  }else{
    X <- X.ini <- as.matrix(X)
    if(intercept)
      X <- cbind(rep(1,k),X)
  }
	
  ##optimize depending on if mu is specified. 
  if(muknown){
    r<-0
    est <- optim(mean(log(V)),function(x){gr.ll.muknown(exp(x),y,V,mu,type)},control=list(fnscale=-1, maxit = 1000),method="BFGS",hessian=T)
    Ahat <- exp(est$par)
    Betahat<-NA
    Betahatvar<-NA
  }else{
    r <- dim(X)[2]
    est <- optim(mean(log(V)), function(x){gr.ll.muunknown(exp(x[1]),y,V,X,type)},control=list(fnscale=-1, maxit = 1000),method="BFGS", hessian=T, gr = function(x){derval(x,y,V,X)})
    Ahat <- exp(est$par[1])
  }

  if(est$convergence == 1)
    warning("Algorithm did not converge")
  
  ##calculate estimates
  est.var <- -1/est$hessian
  ninfo <- -1*est$hessian
  Avar<-est.var*Ahat^2
  Bhat<-V/(V+Ahat)
  a0<-ninfo/Bhat
  a1<-ninfo/(1-Bhat)
  v<-Bhat*(1-Bhat)/(a1+a0+1)
  seB<-sqrt(v)
  skewB<-2*(a0-a1)*sqrt(1+a0+a1)/(sqrt(a0*a1)*(2+a0+a1))
  LCLB <- qbeta((1-confidence.lvl)/2,a1,a0)
  UCLB <- qbeta(1 - (1-confidence.lvl)/2,a1,a0)
  
  ##if mu is estimated use different formula for posterior sd and estimate mu
  if(muknown){
    shat<-sqrt(V*(1-Bhat)+v*(y-mu)^2)
  }else{
    DVAhat<-diag(V+Ahat)
    DVAhati <- diag(1/(V+Ahat))
    Betahatvar <- chol2inv(chol(t(X)%*%DVAhati%*%X))
    BetahatvarmX <- Betahatvar%*%t(X)
    Betahat <- BetahatvarmX%*%DVAhati%*%y
    BetahatSE <- sqrt(diag(Betahatvar))
    mu<-X%*%Betahat
    Z <- sqrt(DVAhati)
    Phat<-Z%*%X%*%BetahatvarmX%*%Z
    p<-diag(Phat)
    shat<-sqrt((1-(1-p)*Bhat)*V+v*(y-mu)^2)
  }
  
  ##calculate some other moments
  thetahat<-(1-Bhat)*y+Bhat*mu
  mu3<-skewB*seB^3
  skew<-((mu-y)^3*mu3-3*V*(mu-y)*v)/shat^3
  sqrtV<-sqrt(V)
  hmB <- length(Bhat)/sum(1/Bhat)

  ## calculate CIs using skewed normal
  ## TODO: use apply not loop

  if(normal.CI){
    snparam <- NULL
    skewedmat <- matrix(nrow = length(thetahat), ncol = 3)
    skewedmat[,1] <- qnorm((1-confidence.lvl)/2,thetahat,shat)
    skewedmat[,2] <- thetahat
    skewedmat[,3] <- qnorm(1-(1-confidence.lvl)/2,thetahat,shat)
  }else{
    snparam <- lapply(1:length(thetahat), function(i){as.numeric(gr.cp.to.dp(c(thetahat[i],shat[i],sign(skew[i])*min(abs(skew[i]),0.94))))})
    tmp <- lapply(1:length(thetahat), function(i){c(qsn((1-confidence.lvl)/2,snparam[[i]][1], snparam[[i]][2], snparam[[i]][3], engine = "biv.nt.prob"),thetahat[i],qsn(1-(1-confidence.lvl)/2,snparam[[i]][1], snparam[[i]][2], snparam[[i]][3], engine = "biv.nt.prob"))})
    skewedmat <- as.matrix(do.call("rbind", tmp))
  }
  
  ## return output
  ## TODO: discuss with Tak and sync output
  output<- list(sample.mean=y,se=se,prior.mean=prior.mean, prior.mean.hat = mu,shrinkage=Bhat,sd.shrinkage=seB,post.intv.low=skewedmat[,1],post.mean=thetahat,post.intv.upp=skewedmat[,3],post.sd=shat,model="gr",X=X.ini,beta.new=Betahat,beta.var=Betahatvar,intercept=intercept,a.new=log(Ahat),a.var=est.var, confidence.lvl = confidence.lvl, weight = NA, snparam = snparam)
  return(output)
}

##log adjusted posterior for known mu
gr.ll.muknown<-function(A,y,V,mu,type){
  lla<-type*log(A)+sum(dnorm(x=y,mean=mu,sd=sqrt(V+A),log=TRUE))
  return(lla)
}

##log adjusted posterior for unknown mu
gr.ll.muunknown<-function(A,y,V,X,type){
  DVAi <- diag(1/(V+A))
  exp1 <- t(X)%*%DVAi%*%X
  BetaA <- try(chol2inv(chol(exp1)) %*% t(X) %*% DVAi %*% y, silent=TRUE)
  if (class(BetaA) == "try-error") lla <- -10^6
  else lla <- type * log(A) + sum(dnorm(x=y, mean=X %*% BetaA, sd=sqrt(V+A), log=TRUE)) - 1/2*log(det(exp1))
  return(lla)
}

##this function is a slightly modified version of cp.to.dp in the sn package
##Author: Adelchi Azzalini
##Webpage: http://azzalini.stat.unipd.it/SN/
gr.cp.to.dp <- function (param) {
  b <- sqrt(2/pi)
  m <- length(param) - 2
  gamma1 <- param[m + 2]
  if (abs(gamma1) > 0.995271746431) 
    gamma1 <- 0.9952  ##this line was altered from original
  A <- sign(gamma1) * (abs(2 * gamma1/(4 - pi)))^(1/3)
  delta <- A/(b * sqrt(1 + A^2))
  lambda <- delta/sqrt(1 - delta^2)
  E.Z <- b * delta
  sd.Z <- sqrt(1 - E.Z^2)
  location <- param[1:m]
  location[1] <- param[1] - param[m + 1] * E.Z/sd.Z
  scale <- param[m + 1]/sd.Z
  dp <- c(location, scale, lambda)
  names(dp)[(m + 1):(m + 2)] <- c("scale", "shape")
  if (m == 1) 
    names(dp)[1] <- "location"
  dp
}


derval <- function(alpha,y,V,X){
  
  A = exp(alpha)
  wv = 1/(V+A)
  Wm = diag(wv)
  tXwm <- t(X)%*%Wm
  Sigmam = chol2inv(chol(tXwm%*%X))
  Betahat <- Sigmam%*%tXwm%*%y
  res <- (y-X%*%Betahat)
  exp1 <- Sigmam%*%t(X)%*%Wm^2
  l2 = log(A) - 1/2*sum(log(V+A)) + 1/2*log(det(Sigmam)) - 1/2*sum(wv*res^2)

  dlralphaBEVAL = 1 - A/2*sum(wv) + A/2*tr(exp1%*%X) + A/2*sum(wv^2*res^2)
  dbetahatA = exp1%*%res
  dl2alpha = dlralphaBEVAL + A*sum(wv*res*X%*%dbetahatA)
  return(dl2alpha=dl2alpha)
}


tr <- function(M){
  sum(diag(M))
}

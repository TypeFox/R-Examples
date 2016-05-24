ttest.Gibbs.AR = function(before,after,iterations=1000,areaNull=c(-.2,.2),leftSided = TRUE,treat=NULL,r.scale=1,alphaTheta=1,betaTheta=5,
                          sdMet=.3, progress=TRUE,return.chains=FALSE, return.onesided = FALSE)
{
    
  y = c(before,after)
  N = length(y)
  
  if(is.null(treat))
  {
    treat = c(rep(-0.5,length(before)),rep(0.5,length(after)))
  }else{
    if(length(treat)!=length(y))
    {
      stop("Invalid condition vector: treat.")
    }
  }
  
  iterations = as.integer(iterations)
  
  if(progress){
    progress = round(iterations/100)
    pb = txtProgressBar(min = 0, max = as.integer(iterations), style = 3) 
  }else{ 
    pb=NULL 
  }
  
  pbFun = function(samps){ if(progress) setTxtProgressBar(pb, samps)}
  
  miss = as.integer(which(is.na(y)) - 1)
  Nmiss = as.integer(length(miss))
  
  yimp = y
  yimp[miss+1] =  mean(y, na.rm =TRUE)
  
  out = .Call("RgibbsTwoSampleAR", as.numeric(yimp), N, treat, r.scale, alphaTheta, betaTheta, 
                 areaNull[1], areaNull[2], iterations, leftSided, sdMet,  miss, Nmiss, progress, 
                 pbFun, new.env(), package="BayesSingleSub")
  
  if(progress) close(pb)
  
  chains = out[[1]]
  postp = out[[2]]
  
  dim(chains) = c(iterations,7)
  chains = data.frame(chains)
  colnames(chains) = c("mu0","beta","ldens","sig2","g","rho","nullArea")
  chains$delta = chains$beta/sqrt(chains$sig2)
  
  logdens = logMeanExpLogs(chains$ldens)
  nulllogdens = dcauchy(0,log=TRUE) - log(r.scale)
  logbf = logdens - nulllogdens
  
  logbfArea = log(mean(chains$nullArea)) - log(1-mean(chains$nullArea)) - 
    (log(diff(pcauchy(areaNull,scale=r.scale))) - log(1 - diff(pcauchy(areaNull,scale=r.scale))))
  
  acc = mean(diff(chains$rho)!=0)
  cat("\n","rho acceptance rate:",acc,"\n")
  
  if(return.onesided == FALSE){
  
    if(return.chains)
    {
      return(list(logbf=logbf,chains=mcmc(chains[,c(1,8,4,5,6)]),acc=acc,logbfArea=logbfArea))
    }else{
      return(logbf)
    }
    
  }
  
  if(return.onesided == TRUE){
    
    logbfOne = log(.5) - postp + logbf
    
    if(return.chains)
    {
      return(list(logbf=logbf,chains=mcmc(chains[,c(1,8,4,5,6)]),acc=acc,logbfArea=logbfArea, logbfOnesided = logbfOne))
    }else{
      return(list(logbf = logbf, logbfOnesided = logbfOne))
    }
    
  }
  
}

ttest.MCGQ.AR = function(before,after,iterations=1000,treat=NULL,method="MC",r.scale=1,alphaTheta=1,betaTheta=5)
{
  
  y = c(before,after)
  
  miss = as.integer(which(is.na(y)))
  Nmiss = as.integer(length(miss))
  
  if(Nmiss > 0){
    y =  y[-miss]
  }
   
  N = length(y)
  
  distMat = abs(outer(1:(N+Nmiss),1:(N+Nmiss),'-')) 
  if(Nmiss > 0){
  distMat = distMat[-miss,-miss]
  }
  oneVec = matrix(1,nrow=N)
  Jn = matrix(1,N,N)
  
  
  if(is.null(treat))
  {
    treat = c(rep(-0.5,length(before)),rep(0.5,length(after)))
  }else{
    if(length(treat)!= (length(y) + Nmiss))
    {
      stop("Invalid condition vector: treat.")
    }
  }
  
  if(Nmiss > 0){
  treat = treat[-miss]
  }

  iterations = as.integer(iterations)
  
  
  vmlike0 = Vectorize(mlike.null.AR,"theta")
  mlike0.gq = log(integrate(vmlike0,lower=0,upper=1,y=as.numeric(y),N=N,tr=treat,oneVec=oneVec,Jn = Jn, alphaTheta=alphaTheta, betaTheta=betaTheta, distMat=distMat)[[1]])
  
    
  if(method=="MC")
  {
    logbf = .Call("RMCTwoSampleAR", as.numeric(y), N, as.integer(distMat), treat, r.scale, alphaTheta, betaTheta, iterations, package="BayesSingleSub") - mlike0.gq
  }else if(method=="GQ")
  {
    vmlike.alt.gtheta.AR = Vectorize(mlike.alt.gtheta.AR,"theta")
    vmlike.alt.g.AR = Vectorize(mlike.alt.g.AR,"g")
    mlike1.gq = log(integrate(vmlike.alt.g.AR,lower=0,upper=Inf,y=y,N=N,tr=treat,oneVec=oneVec,Jn = Jn,alphaTheta=alphaTheta,betaTheta=betaTheta,fun1=vmlike.alt.gtheta.AR,distMat=distMat)[[1]])
    logbf = mlike1.gq - mlike0.gq
  }else{
    stop(paste("Invalid method specified:",method))
  }
  
  return(-logbf) 
  
 }


mlike.null.AR = function(theta,y,tr,N=length(y),oneVec=matrix(y*0+1,ncol=1),Jn = oneVec%*%t(oneVec),alphaTheta=1,betaTheta=5,...)
{
  psifunc = function(theta,distMat) theta^distMat/(1-theta^2)
  
  invpsi = solve(psifunc(theta,...))
  invpsi0 = invpsi - invpsi %*% Jn %*% invpsi/as.vector(t(oneVec)%*%invpsi%*%oneVec)
  ret = - 0.5 * log(t(oneVec)%*%invpsi%*%oneVec) +
    -(N-1)/2 * log(t(y)%*%invpsi0%*%y) + 0.5 * determinant(invpsi,logarithm=TRUE)$modulus + dbeta(theta,alphaTheta,betaTheta,log=TRUE)
  exp(ret)
}

mlike.alt.gtheta.AR = function(theta,g,y,tr,N=length(y),oneVec=matrix(y*0+1,ncol=1),Jn = oneVec%*%t(oneVec),alphaTheta=1,betaTheta=5,rscale=1,...)
{
  psifunc = function(theta,distMat) theta^distMat/(1-theta^2)
  invpsi = solve(psifunc(theta,...))
  opsio = as.vector(t(oneVec)%*%invpsi%*%oneVec)
  invpsi0 = invpsi - invpsi%*%Jn%*%invpsi/opsio 
  tpsi0t = as.vector(t(tr)%*%invpsi0%*%tr) + 1/g
  invpsi1 = invpsi0 - invpsi0%*%tr%*%t(tr)%*%invpsi0/tpsi0t  
  devs = t(y)%*%invpsi1%*%y
  exp(
    -((N-1)/2)*log(devs) + 
      (-0.5)*log(opsio) + (-0.5)*log(tpsi0t) +
      0.5*determinant(invpsi,logarithm=TRUE)$modulus +
      -0.5*log(g) + dbeta(theta,alphaTheta,betaTheta,log=TRUE) + log(dinvgamma(g,0.5,rscale^2/2))
  )
}

mlike.alt.g.AR = function(g,y,tr,N=length(y),psifunc,oneVec=matrix(y*0+1,ncol=1),Jn = oneVec%*%t(oneVec),alphaTheta=1,betaTheta=5,rscale=1,fun1,...)
{
  integrate(fun1,lower=0,upper=1,g=g,y=y,N=N,tr=tr,oneVec=oneVec,Jn = oneVec%*%t(oneVec),rscale=rscale,alphaTheta=alphaTheta,betaTheta=betaTheta,...)[[1]]
}







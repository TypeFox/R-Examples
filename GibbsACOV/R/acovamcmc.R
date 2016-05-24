#' Gibbs sampler for one-way mixed-effects ANCOVA models using flat priors.
#' @param Y Vector of reponses of n subjects
#' @param trt Vector of categorical factor levels for n subjects
#' @param X Design matrix with dimension (n x p) where p is the number of continuous predictors (for ANOVA, p = 1 for include grand mean)
#' @param nochn Number of chains to test convergence of the Gibbs sampler
#' @param numIter Number of iterations in the Gibbs sampler
#' @param initval Matrix of initial values for Gibbs sampler with dimension (nochn, (p + nlevels(trt) + 2))
#' @param credint Coverage probability for parameter credible intervals
#' @param Rthresh Gelman-Rubin diagnostic for test of convergence
#' @return S3 \code{acovamcmc} object; a list consisting of
#'   \item{beta}{values of regression coefficients for each iteration}
#'   \item{sig2a}{values of mixed-effect variance for each iteration}
#'   \item{sig2e}{values of error variance for each iteration}
#'   \item{Credible_Interval}{lower bound, point estimate, and upper bound for parameters}
#'   \item{Credible_Interval_Coverage}{coverage percentage for credible intervals}
#'   \item{Convergence_Diag}{status of Gibbs sampler convergence using threshold set for Gelman and Rubin's diagnostic}
#'   \item{Gelman_Rubin_Threshold}{threshold set for Gelman and Rubin's diagnostic}
#'   \item{Iterations}{number of iterations of Gibbs sampler}
#'   \item{Run_Time}{total elapsed seconds}
#' @examples
#'   library(GibbsACOV)
#'   # ANCOVA with 2 continuous predictors and 5 factor levels
#'   init1 <- c(rep(0,7), 1, 1)
#'   init2 <- c(rnorm(7), rgamma(2,2,1))
#'   init3 <- c(rnorm(7), rgamma(2,2,1))
#'   init4 <- c(rnorm(7), rgamma(2,2,1))
#'   initval <- rbind(init1, init2, init3, init4)
#'   acovamcmc(corn$yield, corn$variety, cbind((corn$nitrogen)^2, corn$nitrogen), 4, 1000 , initval)
#'   # ANOVA with grand mean parameterization and 12 factor levels
#'   init1 <- c(rep(0,13), 1, 1)
#'   init2 <- c(rnorm(13), rgamma(2,2,1))
#'   init3 <- c(rnorm(13), rgamma(2,2,1))
#'   init4 <- c(rnorm(13), rgamma(2,2,1))
#'   initval <- rbind(init1, init2, init3, init4)
#' acovamcmc(cs$rate, factor(cs$hospital), matrix(1,length(cs$hospital),1), 4, 1000, initval)
#' @references Gelman, A and Rubin, DB (1992) Inference from iterative simulation using multiple sequences, Statistical Science, 7, 457-511.

acovamcmc <- function(Y,trt,X,nochn,numIter,initval,credint=0.95,Rthresh=1.1){
  ptm=proc.time() 
  if(length(Y)!=length(trt)||length(Y)!=nrow(X)) stop("data lengths differ")
  ftrt=factor(trt)
  #ntrt is the number of treatment
  ltrt=levels(ftrt)
  ntrt=length(ltrt)
  p=ncol(X)
  numpar=p+ntrt+2
  numBurn = numIter/2
  if(nrow(initval)!=nochn||ncol(initval)!=numpar) stop("initial value dimension error")
  nGlob=length(Y)
  beta=array(1,c(nochn,numIter,p))
  alpha=array(1,c(nochn,numIter,ntrt))
  sig2a=matrix(1,nrow=nochn,ncol=numIter)
  sig2e=matrix(1,nrow=nochn,ncol=numIter)
  #set initial value, initval is the initial value
  for(i in 1:nochn){
    beta[i,1,]=initval[i,1:p]
    alpha[i,1,]=initval[i,(p+1):(p+ntrt)]
    sig2a[i,1]=initval[i,p+ntrt+1]
    sig2e[i,1]=initval[i,p+ntrt+2]
  }
  n=numeric()
  for(j in 1:ntrt)
    n[j]=sum(trt==ltrt[j])
  #the following is the "Gibbs sampler"
  for(chn in 1:nochn){
  for(i in 2:numIter){
    #first sample beta
    mu=ginv(t(X)%*%X)%*%t(X)%*%(Y-alpha[chn,i-1,trt])
    Sigma=ginv(t(X)%*%X/sig2e[chn,i-1])
    beta[chn,i,]=mvrnorm(1,mu,Sigma)
    #then sample alpha
    for(j in 1:ntrt){
      sig2=1/(1/sig2a[chn,i-1]+n[j]/sig2e[chn,i-1])
      if(dim(beta)[3]==1)  mu=sig2*sum(Y[trt==ltrt[j]]-X[trt==ltrt[j],]*beta[chn,i,])/sig2e[chn,i-1]
      else mu=sig2*sum(Y[trt==ltrt[j]]-X[trt==ltrt[j],]%*%beta[chn,i,])/sig2e[chn,i-1]
      alpha[chn,i,j]=rnorm(1,mu,sd=sqrt(sig2))
    }
    #then sample sig2a
    sig2a[chn,i]=1/rgamma(1,shape=1+ntrt/2,scale=2/sum(alpha[chn,i,]^2))
    #then sample sig2e
    if(dim(beta)[3]==1)
    sig2e[chn,i]=1/rgamma(1,shape=nGlob/2+1,scale=2/sum((Y-X*beta[chn,i,]-alpha[chn,i,trt])^2))
    else
    sig2e[chn,i]=1/rgamma(1,shape=nGlob/2+1,scale=2/sum((Y-X%*%beta[chn,i,]-alpha[chn,i,trt])^2))
    }
  }
  credlevel=1-credint
  if(ncol(X)>1){
    beta_hat=apply(beta[1,(numBurn+1):(numIter),],2,mean)
    beta_ll=apply(beta[1,(numBurn+1):(numIter),],2,quantile,credlevel/2)
    beta_uu=apply(beta[1,(numBurn+1):(numIter),],2,quantile,1-credlevel/2)
  }
  if(ncol(X)==1){
    beta_hat=mean(beta[1,(numBurn+1):(numIter),])
    beta_ll=quantile(beta[1,(numBurn+1):(numIter),],credlevel/2)
    beta_uu=quantile(beta[1,(numBurn+1):(numIter),],1-credlevel/2)
  }
  beta_CI=cbind(beta_ll,beta_hat,beta_uu)
  sig2a_hat=mean(sig2a[1,(numBurn+1):(numIter)])
  sig2a_ll=quantile(sig2a[1,(numBurn+1):(numIter)],credlevel/2)
  sig2a_uu=quantile(sig2a[1,(numBurn+1):(numIter)],1-credlevel/2)
  sig2e_hat=mean(sig2e[1,(numBurn+1):(numIter)])
  sig2e_ll=quantile(sig2e[1,(numBurn+1):(numIter)],credlevel/2)
  sig2e_uu=quantile(sig2e[1,(numBurn+1):(numIter)],1-credlevel/2)
  CI=rbind(beta_CI,c(sig2a_ll,sig2a_hat,sig2a_uu),c(sig2e_ll,sig2e_hat,sig2e_uu))
  colnames(CI)=c("lower","est","upper")
  nbetas=ncol(X)
  namesCI=list(NA,nbetas+2)
  for(l in 1:nbetas){
    namesCI[[l]]=paste("beta",l)
  }
  namesCI[[nbetas+1]]="sigma_a^2"
  namesCI[[nbetas+2]]="sigma_e^2"
  rownames(CI)=namesCI
 
  #Gelman-Rubin diag
  cnvgcnt=0
  numBurn=numIter/2
  #test convergence of beta
  s=matrix(0,nrow=p,ncol=nochn)
  b=matrix(0,nrow=p,ncol=nochn)
  B=numeric()
  w=numeric()
  v=numeric()
  R=numeric()
  for(i in 1:p){
    for(j in 1:nochn){
      s[i,j]=var(beta[j,(numBurn+1):numIter,i])
      b[i,j]=(mean(beta[j,(numBurn+1):numIter,i])-mean(beta[,(numBurn+1):numIter,i]))^2
    }
    B[i]=sum(b[i,])*numBurn/(nochn-1)
    w[i]=sum(s[i,])/nochn
    v[i]=(1-1/numBurn)*w[i]+B[i]/numBurn
    R[i]=sqrt(v[i]/w[i])
  }
  cnvgcnt=cnvgcnt+sum(R<Rthresh)
  #test convergence of alpha
  s=matrix(0,nrow=ntrt,ncol=nochn)
  b=matrix(0,nrow=ntrt,ncol=nochn)
  B=numeric()
  w=numeric()
  v=numeric()
  R=numeric()
  for(i in 1:ntrt){
    for(j in 1:nochn){
      s[i,j]=var(alpha[j,(numBurn+1):numIter,i])
      b[i,j]=(mean(alpha[j,(numBurn+1):numIter,i])-mean(alpha[,(numBurn+1):numIter,i]))^2
    }
    B[i]=sum(b[i,])*numBurn/(nochn-1)
    w[i]=sum(s[i,])/nochn
    v[i]=(1-1/numBurn)*w[i]+B[i]/numBurn
    R[i]=sqrt(v[i]/w[i])
  }
  cnvgcnt=cnvgcnt+sum(R<Rthresh)  
  #test convergence of sig2a
  s=matrix(0,nrow=1,ncol=nochn)
  b=matrix(0,nrow=1,ncol=nochn)
  B=numeric()
  w=numeric()
  v=numeric()
  R=numeric()
  for(i in 1:1){
    for(j in 1:nochn){
      s[i,j]=var(sig2a[j,(numBurn+1):numIter])
      b[i,j]=(mean(sig2a[j,(numBurn+1):numIter])-mean(sig2a[,(numBurn+1):numIter]))^2
    }
    B[i]=sum(b[i,])*numBurn/(nochn-1)
    w[i]=sum(s[i,])/nochn
    v[i]=(1-1/numBurn)*w[i]+B[i]/numBurn
    R[i]=sqrt(v[i]/w[i])
  }
  cnvgcnt=cnvgcnt+sum(R[i]<Rthresh)  
  #test convergence of sig2e
  s=matrix(0,nrow=1,ncol=nochn)
  b=matrix(0,nrow=1,ncol=nochn)
  B=numeric()
  w=numeric()
  v=numeric()
  R=numeric()
  for(i in 1:1){
    for(j in 1:nochn){
      s[i,j]=var(sig2e[j,(numBurn+1):numIter])
      b[i,j]=(mean(sig2e[j,(numBurn+1):numIter])-mean(sig2e[,(numBurn+1):numIter]))^2
    }
    B[i]=sum(b[i,])*numBurn/(nochn-1)
    w[i]=sum(s[i,])/nochn
    v[i]=(1-1/numBurn)*w[i]+B[i]/numBurn
    R[i]=sqrt(v[i]/w[i])
  }
  cnvgcnt=cnvgcnt+sum(R[i]<Rthresh)  
  if(cnvgcnt>=numpar){cvg=paste("Convergent after",numIter,"iterations.")
  	} else cvg=paste("Not convergent after",numIter,"iterations.")

  runtime=proc.time()-ptm

  result = list(beta=beta[1,,],sig2a=sig2a[1,],sig2e=sig2e[1,],Credible_Interval=CI,Credible_Interval_Coverage=credint,Convergence_Diag=cvg,Gelman_Rubin_Threshold=Rthresh,Iterations=numIter,Run_Time=runtime)
  class(result)="acovamcmc"
  return(result)
}

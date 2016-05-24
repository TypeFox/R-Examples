drvts <-
function(bug, sims, ymean, hmean=NULL, iter=NULL){
  if(is.null(colnames(sims)))
    stop("columns of sims can not be NULL. names should correspond to parameters")
  if(class(bug)!="tsbugs")
    stop("bug must be a object of class tsbugs")
  
  theta<-theta.it(bug,sims)
  k<-nrow(sims)
  stoc<-NULL
  dd<-nodes(bug, part="prior")
  dd<-subset(dd, stoc==1)
  
  phi.prior<-theta$phi
  phi.prior[]<-0
  j<-match("phi1",dd$name)
  if(!is.na(j)){
    f<-match.fun(dd$dist[j])
    phi.prior<-f(theta$phi,dd$param1[j],dd$param2[j],log=TRUE)
    phi.prior[theta$phi==0]<-0              #set densitites to zero for ones not in model
  }
  j<-match("phi0",dd$name)
  if(!is.na(j)){
    f<-match.fun(dd$dist[j])
    phi.prior[,"phi0"]<-f(theta$phi[,"phi0"],dd$param1[j],dd$param2[j],log=TRUE)
  }
  phi.prior<-apply(phi.prior,1,sum) #sum for each simulation
  
  j<-match("isig02",dd$name)
  f<-match.fun(dd$dist[j])
  isig02.prior<-f(theta$isig02,dd$param1[j],dd$param2[j], log=TRUE)
  
  j<-match("ilambda2",dd$name)
  f<-match.fun(dd$dist[j])
  ilambda2.prior<-f(theta$ilambda2,dd$param1[j],dd$param2[j], log=TRUE)
  
  j<-match("epsilon",dd$name)
  f<-match.fun(dd$dist[j])
  epsilon.prior<-f(theta$epsilon, dd$param1[j],dd$param2[j], log=TRUE)
  
  beta.prior<-dnorm(theta$beta, 0, sqrt(1/theta$ilambda2), log=TRUE)
  beta.prior<-apply(beta.prior,1,sum)
  delta.prior<-dbinom(theta$delta, 1, theta$epsilon, log=TRUE)
  delta.prior<-apply(delta.prior,1,sum)
  
  #likelihood
  mod.lik<-matrix(NA,k,1)
  for(i in 1:k){
    mod.lik[i,]<-tslogl(bug, ymean=ymean[i,], sigma=sqrt(exp(hmean[i,])))
  }
  
  #posterior
  post<-phi.prior+isig02.prior+beta.prior+delta.prior+epsilon.prior+ilambda2.prior+mod.lik
  
  if(is.numeric(iter)){
    disp<-list(phi=phi.prior[iter], 
               epsilon=epsilon.prior[iter],
               isig02=isig02.prior[iter], 
               beta=beta.prior[iter], 
               delta=delta.prior[iter], 
               ilambda2=ilambda2.prior[iter],
               lik=mod.lik[iter], post=post[iter])
    print(disp)
  }
  if(is.null(iter)){
    return(post)
  }
}

dcvts <-
function(bug, sims, ymean, hmean=NULL, iter=NULL){
  if(is.null(colnames(sims)))
    stop("columns of sims can not be NULL. names should correspond to parameters")
  if(class(bug)!="tsbugs")
    stop("bug must be a object of class tsbugs")
  
  #components of theta
  theta<-theta.it(bug,sims)
  k<-nrow(sims)
  stoc<-NULL
  dd<-nodes(bug, part="prior")
  dd<-subset(dd, stoc==1)
  #m<-bug$info$mlag  
  #nodes(bug,"prior")
  
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
  
  sigma.prior<-rep(0,k)
  sigma2.prior<-rep(0,k)
  isigma2.prior<-rep(0,k)
  if(!is.null(theta$sigma)){
    j<-match("sigma",dd$name)
    f<-match.fun(dd$dist[j])
    sigma.prior<-f(theta$sigma,dd$param1[j],dd$param2[j], log=TRUE)
  }
  if(!is.null(theta$sigma2)){
    j<-match("sigma2",dd$name)
    f<-match.fun(dd$dist[j])
    sigma2.prior<-f(theta$sigma2,dd$param1[j],dd$param2[j], log=TRUE)
  }
  if(!is.null(theta$isigma2)){
    j<-match("isigma2",dd$name)
    f<-match.fun(dd$dist[j])
    isigma2.prior<-f(theta$isigma2,dd$param1[j],dd$param2[j], log=TRUE)
  }
  
  if(is.null(theta$sigma) & !is.null(theta$sigma2)){
    theta$sigma<-sqrt(theta$sigma2)
  }
  if(is.null(theta$sigma) & !is.null(theta$isigma2)){
    theta$sigma<-1/sqrt(theta$isigma2)
  }
  #likelihood
  mod.lik<-matrix(NA,k,1)
  for(i in 1:k){
    mod.lik[i,]<-tslogl(bug, ymean=ymean[i,], sigma=theta$sigma[i])
  }
  
  #posterior
  post<-phi.prior+sigma.prior+sigma2.prior+isigma2.prior+mod.lik
  
  if(is.numeric(iter)){
    disp<-list(phi=phi.prior[iter], 
               sigma=sigma.prior[iter], 
               sigma2=sigma2.prior[iter], 
               isigma2=isigma2.prior[iter], 
               lik=mod.lik[iter], post=post[iter])
    print(disp)
  }
  if(is.null(iter)){
    return(post)
  }
}

dsvts <-
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
  
  psi.prior<-theta$psi
  psi.prior[]<-0
  j<-match("psi1",dd$name)
  if(!is.na(j)){
    f<-match.fun(dd$dist[j])
    psi.prior<-f(theta$psi,dd$param1[j],dd$param2[j],log=TRUE)
    psi.prior[theta$psi==0]<-0              #set densitites to zero for ones not in model
  }
  j<-match("psi0",dd$name)
  if(!is.na(j)){
    f<-match.fun(dd$dist[j])
    psi.prior[,"psi0"]<-f(theta$psi[,"psi0"],dd$param1[j],dd$param2[j],log=TRUE)
  }
  psi.prior<-apply(psi.prior,1,sum) #sum for each simulation
  
  psi.star.prior<-theta$psi.star
  psi.star.prior[]<-0
  j<-match("psi1.star",dd$name)
  if(!is.na(j)){
    f<-match.fun(dd$dist[j])
    psi.star.prior<-f(theta$psi.star,dd$param1[j],dd$param2[j],log=TRUE)
    psi.star.prior[theta$psi.star==0]<-0              #set densitites to zero for ones not in model
  }
  j<-match("psi0.star",dd$name)
  if(!is.na(j)){
    f<-match.fun(dd$dist[j])
    psi.star.prior[,"psi0.star"]<-f(theta$psi.star[,"psi0.star"],dd$param1[j],dd$param2[j],log=TRUE)
  }
  psi.star.prior<-apply(psi.star.prior,1,sum) #sum for each simulation
  
  j<-match("itau2",dd$name)
  f<-match.fun(dd$dist[j])
  itau2.prior<-f(theta$itau2,dd$param1[j],dd$param2[j], log=TRUE)
  
  eta.prior<-dnorm(theta$h-hmean, 0, 1/sqrt(theta$itau2), log=TRUE)
  eta.prior<-apply(eta.prior,1,sum)
  
  #likelihood
  mod.lik<-matrix(NA,k,1)
  for(i in 1:k){
     mod.lik[i,]<-tslogl(bug, ymean=ymean[i,], sigma=sqrt(exp(theta$h[i,])))
#      mod.lik[i,]<-tslogl(bug, ymean=ymean[i,], sigma=sqrt(exp(hmean[i,])))
  }

  #posterior
  post<-phi.prior+psi.prior+psi.star.prior+itau2.prior+eta.prior+mod.lik
  
  if(is.numeric(iter)){
    disp<-list(phi=phi.prior[iter], 
               psi=psi.prior[iter], 
               psi.star=psi.star.prior[iter], 
               itau2=itau2.prior[iter], 
               eta=eta.prior[iter], 
               lik=mod.lik[iter], post=post[iter])
    print(disp)
  }
  if(is.null(iter)){
    return(post)
  }
}

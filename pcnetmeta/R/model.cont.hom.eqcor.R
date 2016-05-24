model.cont.hom.eqcor <- function(prior.type="unif",rank.prob=TRUE){
if(prior.type=="unif" & rank.prob){
modelstring<-"
model{
 for(i in 1:len){
  mean[i]~dnorm(theta[i],n[i]/pow(sd[i],2))
  theta[i]<-mu[t[i]]+vi[s[i],t[i]]
 }
 for(j in 1:nstudy){
  vi[j,1:ntrt]~dmnorm(zeros[1:ntrt],T[1:ntrt,1:ntrt])
 }
 for(j in 1:ntrt){
  mu[j]~dnorm(0,0.001)
 }
 for(i in 1:ntrt){
  for(j in 1:ntrt){
   diff[i,j]<-mu[i]-mu[j]
  }
 }
 for(j in 1:ntrt){
  for(k in 1:ntrt){ 
   T[j,k]<-1/sigma^2*ifelse(j==k,diag,offdiag)
  }
 }
 diag<-(1+(ntrt-2)*rho)/(1+(ntrt-2)*rho-(ntrt-1)*rho^2)
 offdiag<-(-rho/(1+(ntrt-2)*rho-(ntrt-1)*rho^2))
 rho~dunif(-1/(ntrt-1),1)
 sigma~dunif(0,c)
 rk[1:ntrt]<-(ntrt+1-rank(mu[]))*ifelse(higher.better,1,0)+(rank(mu[]))*ifelse(higher.better,0,1)
 for(i in 1:ntrt){
  rank.prob[1:ntrt,i]<-equals(rk[],i)
 }
}
"
}

if(prior.type=="unif" & !rank.prob){
modelstring<-"
model{
 for(i in 1:len){
  mean[i]~dnorm(theta[i],n[i]/pow(sd[i],2))
  theta[i]<-mu[t[i]]+vi[s[i],t[i]]
 }
 for(j in 1:nstudy){
  vi[j,1:ntrt]~dmnorm(zeros[1:ntrt],T[1:ntrt,1:ntrt])
 }
 for(j in 1:ntrt){
  mu[j]~dnorm(0,0.001)
 }
 for(i in 1:ntrt){
  for(j in 1:ntrt){
   diff[i,j]<-mu[i]-mu[j]
  }
 }
 for(j in 1:ntrt){
  for(k in 1:ntrt){ 
   T[j,k]<-1/sigma^2*ifelse(j==k,diag,offdiag)
  }
 }
 diag<-(1+(ntrt-2)*rho)/(1+(ntrt-2)*rho-(ntrt-1)*rho^2)
 offdiag<-(-rho/(1+(ntrt-2)*rho-(ntrt-1)*rho^2))
 rho~dunif(-1/(ntrt-1),1)
 sigma~dunif(0,c)
}
"
}

if(prior.type=="invgamma" & rank.prob){
modelstring<-"
model{
 for(i in 1:len){
  mean[i]~dnorm(theta[i],n[i]/pow(sd[i],2))
  theta[i]<-mu[t[i]]+vi[s[i],t[i]]
 }
 for(j in 1:nstudy){
  vi[j,1:ntrt]~dmnorm(zeros[1:ntrt],T[1:ntrt,1:ntrt])
 }
 for(j in 1:ntrt){
  mu[j]~dnorm(0,0.001)
 }
 for(i in 1:ntrt){
  for(j in 1:ntrt){
   diff[i,j]<-mu[i]-mu[j]
  }
 }
 for(j in 1:ntrt){
  for(k in 1:ntrt){ 
   T[j,k]<-1/sigma^2*ifelse(j==k,diag,offdiag)
  }
 }
 diag<-(1+(ntrt-2)*rho)/(1+(ntrt-2)*rho-(ntrt-1)*rho^2)
 offdiag<-(-rho/(1+(ntrt-2)*rho-(ntrt-1)*rho^2))
 rho~dunif(-1/(ntrt-1),1)
 sigma<-1/sqrt(inv.sig.sq)
 inv.sig.sq~dgamma(a,b)
 rk[1:ntrt]<-(ntrt+1-rank(mu[]))*ifelse(higher.better,1,0)+(rank(mu[]))*ifelse(higher.better,0,1)
 for(i in 1:ntrt){
  rank.prob[1:ntrt,i]<-equals(rk[],i)
 }
}
"
}

if(prior.type=="invgamma" & !rank.prob){
modelstring<-"
model{
 for(i in 1:len){
  mean[i]~dnorm(theta[i],n[i]/pow(sd[i],2))
  theta[i]<-mu[t[i]]+vi[s[i],t[i]]
 }
 for(j in 1:nstudy){
  vi[j,1:ntrt]~dmnorm(zeros[1:ntrt],T[1:ntrt,1:ntrt])
 }
 for(j in 1:ntrt){
  mu[j]~dnorm(0,0.001)
 }
 for(i in 1:ntrt){
  for(j in 1:ntrt){
   diff[i,j]<-mu[i]-mu[j]
  }
 }
 for(j in 1:ntrt){
  for(k in 1:ntrt){ 
   T[j,k]<-1/sigma^2*ifelse(j==k,diag,offdiag)
  }
 }
 diag<-(1+(ntrt-2)*rho)/(1+(ntrt-2)*rho-(ntrt-1)*rho^2)
 offdiag<-(-rho/(1+(ntrt-2)*rho-(ntrt-1)*rho^2))
 rho~dunif(-1/(ntrt-1),1)
 sigma<-1/sqrt(inv.sig.sq)
 inv.sig.sq~dgamma(a,b)
}
"
}

if(!is.element(prior.type,c("unif","invgamma"))){
  stop("specified prior type is wrong.")
}

return(modelstring)
}
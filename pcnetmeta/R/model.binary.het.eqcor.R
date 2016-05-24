model.binary.het.eqcor <- function(prior.type="unif",rank.prob=TRUE){
if(prior.type=="unif" & rank.prob){
modelstring<-"
model{
 for(i in 1:len){
  p[i]<-phi(mu[t[i]]+vi[s[i],t[i]])
  r[i]~dbin(p[i],totaln[i])
 }
 for(j in 1:nstudy){
  vi[j,1:ntrt]~dmnorm(zeros[1:ntrt],T[1:ntrt,1:ntrt])
 }
 for(j in 1:ntrt){
  AR[j]<-phi(mu[j]/sqrt(1+pow(sigma[j],2)))
  mu[j]~dnorm(0,0.001)
  sigma[j]~dunif(0,c)
 }
 for(j in 1:ntrt){
  for(k in 1:ntrt){ 
   T[j,k]<-1/sigma[j]*1/sigma[k]*ifelse(j==k,diag,offdiag)
  }
 }
 diag<-(1+(ntrt-2)*rho)/(1+(ntrt-2)*rho-(ntrt-1)*rho^2)
 offdiag<-(-rho/(1+(ntrt-2)*rho-(ntrt-1)*rho^2))
 rho~dunif(-1/(ntrt-1),1)
 for(j in 1:ntrt){
  for(k in 1:ntrt){
   LRR[j,k]<-log(RR[j,k])
   LOR[j,k]<-log(OR[j,k])
   RR[j,k]<-AR[j]/AR[k]
   RD[j,k]<-AR[j]-AR[k]
   OR[j,k]<-AR[j]/(1-AR[j])/AR[k]*(1-AR[k])
  }
 }
 rk[1:ntrt]<-(ntrt+1-rank(AR[]))*ifelse(higher.better,1,0)+(rank(AR[]))*ifelse(higher.better,0,1)
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
  p[i]<-phi(mu[t[i]]+vi[s[i],t[i]])
  r[i]~dbin(p[i],totaln[i])
 }
 for(j in 1:nstudy){
  vi[j,1:ntrt]~dmnorm(zeros[1:ntrt],T[1:ntrt,1:ntrt])
 }
 for(j in 1:ntrt){
  AR[j]<-phi(mu[j]/sqrt(1+pow(sigma[j],2)))
  mu[j]~dnorm(0,0.001)
  sigma[j]~dunif(0,c)
 }
 for(j in 1:ntrt){
  for(k in 1:ntrt){ 
   T[j,k]<-1/sigma[j]*1/sigma[k]*ifelse(j==k,diag,offdiag)
  }
 }
 diag<-(1+(ntrt-2)*rho)/(1+(ntrt-2)*rho-(ntrt-1)*rho^2)
 offdiag<-(-rho/(1+(ntrt-2)*rho-(ntrt-1)*rho^2))
 rho~dunif(-1/(ntrt-1),1)
 for(j in 1:ntrt){
  for(k in 1:ntrt){
   LRR[j,k]<-log(RR[j,k])
   LOR[j,k]<-log(OR[j,k])
   RR[j,k]<-AR[j]/AR[k]
   RD[j,k]<-AR[j]-AR[k]
   OR[j,k]<-AR[j]/(1-AR[j])/AR[k]*(1-AR[k])
  }
 }
}
"
}

if(prior.type=="invgamma" & rank.prob){
modelstring<-"
model{
 for(i in 1:len){
  p[i]<-phi(mu[t[i]]+vi[s[i],t[i]])
  r[i]~dbin(p[i],totaln[i])
 }
 for(j in 1:nstudy){
  vi[j,1:ntrt]~dmnorm(zeros[1:ntrt],T[1:ntrt,1:ntrt])
 }
 for(j in 1:ntrt){
  AR[j]<-phi(mu[j]/sqrt(1+pow(sigma[j],2)))
  mu[j]~dnorm(0,0.001)
  sigma[j]<-1/sqrt(inv.sig.sq[j])
  inv.sig.sq[j]~dgamma(a,b)
 }
 for(j in 1:ntrt){
  for(k in 1:ntrt){ 
   T[j,k]<-1/sigma[j]*1/sigma[k]*ifelse(j==k,diag,offdiag)
  }
 }
 diag<-(1+(ntrt-2)*rho)/(1+(ntrt-2)*rho-(ntrt-1)*rho^2)
 offdiag<-(-rho/(1+(ntrt-2)*rho-(ntrt-1)*rho^2))
 rho~dunif(-1/(ntrt-1),1)
 for(j in 1:ntrt){
  for(k in 1:ntrt){
   LRR[j,k]<-log(RR[j,k])
   LOR[j,k]<-log(OR[j,k])
   RR[j,k]<-AR[j]/AR[k]
   RD[j,k]<-AR[j]-AR[k]
   OR[j,k]<-AR[j]/(1-AR[j])/AR[k]*(1-AR[k])
  }
 }
 rk[1:ntrt]<-(ntrt+1-rank(AR[]))*ifelse(higher.better,1,0)+(rank(AR[]))*ifelse(higher.better,0,1)
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
  p[i]<-phi(mu[t[i]]+vi[s[i],t[i]])
  r[i]~dbin(p[i],totaln[i])
 }
 for(j in 1:nstudy){
  vi[j,1:ntrt]~dmnorm(zeros[1:ntrt],T[1:ntrt,1:ntrt])
 }
 for(j in 1:ntrt){
  AR[j]<-phi(mu[j]/sqrt(1+pow(sigma[j],2)))
  mu[j]~dnorm(0,0.001)
  sigma[j]<-1/sqrt(inv.sig.sq[j])
  inv.sig.sq[j]~dgamma(a,b)
 }
 for(j in 1:ntrt){
  for(k in 1:ntrt){ 
   T[j,k]<-1/sigma[j]*1/sigma[k]*ifelse(j==k,diag,offdiag)
  }
 }
 diag<-(1+(ntrt-2)*rho)/(1+(ntrt-2)*rho-(ntrt-1)*rho^2)
 offdiag<-(-rho/(1+(ntrt-2)*rho-(ntrt-1)*rho^2))
 rho~dunif(-1/(ntrt-1),1)
 for(j in 1:ntrt){
  for(k in 1:ntrt){
   LRR[j,k]<-log(RR[j,k])
   LOR[j,k]<-log(OR[j,k])
   RR[j,k]<-AR[j]/AR[k]
   RD[j,k]<-AR[j]-AR[k]
   OR[j,k]<-AR[j]/(1-AR[j])/AR[k]*(1-AR[k])
  }
 }
}
"
}

if(!is.element(prior.type,c("unif","invgamma"))){
  stop("specified prior type is wrong.")
}

return(modelstring)
}
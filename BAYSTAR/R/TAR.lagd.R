TAR.lagd<-function(ay,p1,p2,ph.1,ph.2,sig.1,sig.2,thres,lagp1,lagp2,constant=1,d0,thresVar){

loglik<-lik<-pr<-NULL

if (!missing(thresVar)){
for (i in 1:d0){         ## Calculate log-likelihood from d=1 to d=d0
loglik[i]<- TAR.lik(ay,p1,p2,ph.1,ph.2,sig.1,sig.2,i,thres,lagp1,lagp2,constant=constant,thresVar)}}
else  {
for (i in 1:d0){         ## Calculate log-likelihood from d=1 to d=d0
loglik[i]<- TAR.lik(ay,p1,p2,ph.1,ph.2,sig.1,sig.2,i,thres,lagp1,lagp2,constant=constant)}}
lik<- (exp(loglik-max(loglik)))*(rev(c(1:d0))/sum(1:d0)) #give weight for lagd
## Comparing the cdf of delay lag with a random probability
## to determine a new delay lag.
lagd<- (sum((cumsum(lik)/sum(lik))<runif(1, min=0, max=1)))+1
return(lagd)
}


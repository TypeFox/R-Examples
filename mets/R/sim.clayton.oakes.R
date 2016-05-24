##' Simulate observations from the Clayton-Oakes copula model with
##' piecewise constant marginals.
##'
##' @title Simulate from the Clayton-Oakes frailty model
##' @param K Number of clusters
##' @param n Number of observations in each cluster
##' @param eta 1/variance
##' @param beta Effect (log hazard ratio) of covariate
##' @param stoptime Stopping time
##' @param left Left truncation
##' @param pairleft pairwise (1) left truncation or individual (0)
##' @param trunc.prob Truncation probability 
##' @author Klaus K. Holst
##' @export
simClaytonOakes <- function(K,n,eta,beta,stoptime,left=0,pairleft=0,trunc.prob=0.5)  ## {{{ 
{
  ## K antal clustre, n=antal i clustre
  x<-array(c(runif(n*K),rep(0,n*K),rbinom(n*K,1,0.5)),dim=c(K,n,3))
  C<-matrix(stoptime,K,n);
  Gam1 <-matrix(rgamma(K,eta)/eta,K,n)
  temp<-eta*log(-log(1-x[,,1])/(eta*Gam1)+1)*exp(-beta*x[,,3])
  x[,,2]<-ifelse(temp<=C,1,0);
  x[,,1]<-pmin(temp,C)
  minstime <- apply(x[,,1],1,min)  
  ud <- as.data.frame(cbind(apply(x,3,t),rep(1:K,each=n)))  

if (left>0) {
     if (pairleft==1) {
     lefttime <- runif(K)*(stoptime-left)
     left <- rbinom(K,1,trunc.prob) ## not trunation times!
     lefttime <- apply(cbind(lefttime*left,3),1,min)
       trunk <- (lefttime > minstime)
       medleft <- rep(trunk,each=n)
     } else {
###       lefttime <- rexp(n*K)*left
       lefttime <- runif(K)*(stoptime-left)
       left <- rbinom(n*K,1,trunc.prob) ## not trunation times!
       lefttime <- apply(cbind(lefttime*left,3),1,min)
       trunk <- (lefttime > ud[,1])
       medleft <- trunk
     }
  } else { lefttime <- trunk <- rep(0,K);}
  if (pairleft==1) ud <- cbind(ud,rep(minstime,each=n),rep(lefttime,each=n),rep(trunk,each=n))
  else ud <- cbind(ud,rep(minstime,each=n),lefttime,trunk)

###  if (left>0) {
###    lefttime <- rexp(K)*left
###    left <- rbinom(K,1,0.5) ## not trunation times!
###    lefttime <- apply(cbind(lefttime*left,3),1,min)
###    trunk <- (lefttime > minstime)
###    medleft <- rep(trunk,each=n)
###  } else { lefttime <- trunk <- rep(0,K);}
###
###  ud <- cbind(ud,rep(minstime,each=n),rep(lefttime,each=n),rep(trunk,each=n))
  names(ud)<-c("time","status","x1","cluster","mintime","lefttime","truncated")
  return(ud)
} ## }}} 

##' Simulate observations from the Clayton-Oakes copula model with
##' Weibull type baseline and Cox marginals.
##'
##' @title Simulate from the Clayton-Oakes frailty model
##' @param K Number of clusters
##' @param n Number of observations in each cluster
##' @param eta 1/variance
##' @param beta Effect (log hazard ratio) of covariate
##' @param stoptime Stopping time
##' @param weiscale weibull scale parameter 
##' @param weishape weibull shape parameter 
##' @param left Left truncation
##' @param pairleft pairwise (1) left truncation or individual (0)
##' @author Klaus K. Holst 
##' @export
simClaytonOakesWei <- function(K,n,eta,beta,stoptime,
	       weiscale=1,weishape=2,left=0,pairleft=0)
{ ## {{{ 
cat(" not quite \n"); 
## K antal clustre, n=antal i clustre
### K=10; n=2; eta=1; beta=0.3; stoptime=3; lam=0.5; 
### weigamma=2; left=0; pairleft=0
 X <- rbinom(n*K,1,0.5)
 C<-rep(stoptime,n*K);
 Gam1 <-rep(rgamma(K,eta),each=n)
 temp <- rexp(K*n)
### temp <- rweibull(n*K,weishape,scale=weiscale)/(exp(X*beta)*Gam1)
 temp<- (eta*log(eta*temp/(eta*Gam1)+1)/(exp(beta*X)*weiscale^weishape))^{1/weishape}
 status<- ifelse(temp<=C,1,0);
 temp <-   pmin(temp,C)
 xt <- matrix(temp,n,K)
 minstime <- apply(xt,2,min)  
 id=rep(1:K,each=n)
  ud <- cbind(temp,status,X,id)
if (left>0) {
     if (pairleft==1) {
     lefttime <- runif(K)*(stoptime-left)
     left <- rbinom(K,1,0.5) ## not trunation times!
     lefttime <- apply(cbind(lefttime*left,3),1,min)
       trunk <- (lefttime > minstime)
       medleft <- rep(trunk,each=n)
     } else {
###       lefttime <- rexp(n*K)*left
       lefttime <- runif(K)*(stoptime-left)
       left <- rbinom(n*K,1,0.5) ## not trunation times!
       lefttime <- apply(cbind(lefttime*left,3),1,min)
       trunk <- (lefttime > ud[,1])
       medleft <- trunk
     }
  } else { lefttime <- trunk <- rep(0,K);}
  if (pairleft==1) ud <- cbind(ud,rep(minstime,each=n),rep(lefttime,each=n),rep(trunk,each=n))
  else ud <- cbind(ud,rep(minstime,each=n),lefttime,trunk)

  colnames(ud)<-c("time","status","x1","cluster","mintime","lefttime","truncated")
  ud <- data.frame(ud)
  return(ud)
} ## }}} 


## sim.clayton <- function(n=100,K=2,eta=0.5,beta,...) {
##     m <- lvm(T~x)
##     rates <- c(0.3,0.5); cuts <- c(0,5)
##     distribution(m,~) <- coxExponential.lvm(rate=rates,timecut=cuts)    
## }


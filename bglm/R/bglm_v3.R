################
require("mvtnorm")
binommh <- function(y,...) UseMethod("binommh")
binommh.default<-function(y,m,X,b=rep(0,dim(X)[2]),B=diag(rep(10000,dim(X)[2])),N=3000,flag=F,...){
  ########utility functions
  yhat<-function(y,X,bta,W_n=1,mu_n){drop(X%*%bta+(y-mu_n)/W_n)}
  btahat<-function(yhat_n,X,W_n=1){drop(solve(t(X)%*%(W_n*X))%*%t(X)%*%(W_n*yhat_n))}
  log_lik<-function(y,m,X,bta){
    eta_y<-X%*%bta
    sum(dbinom(y,m,drop(exp(eta_y)/(1+exp(eta_y))),log=T)) 
  }
  ###################
  chain<-matrix(rep(NA,dim(X)[2]*N),nrow=N);colnames(chain)<-colnames(X)
  Dev<-rep(NA,N)
  accept <- 0
  ###frcuent cuantities
  Binv<-solve(B)
  ####initials
  ytemp<-y
  ytemp[ytemp==0]<-.5
  W_0<-ytemp*(m-ytemp)/m
  den<-m-ytemp
  den[den==0]<-0.5
  yhat_0<-log(ytemp/(den))
  chain[1,]<-btahat(yhat_0,X,W_0)
  Dev[1]<- -2*log_lik(y,m,X,chain[1,])
  ###proposal kernel
  for(i in 2:N){
    mu_n<-drop(m*exp(X%*%chain[i-1,])/(1+exp(X%*%chain[i-1,])))
    W_n<-mu_n*(m-mu_n)/m
    yhat_n<-yhat(y,X,chain[i-1,],W_n,mu_n)
    
    B.<-solve(Binv+t(X)%*%(W_n*X))
    b.<-drop(B.%*%(Binv%*%b+t(X)%*%(W_n*yhat_n)))
    betaprop<-drop(rmvnorm(1,b.,B.))
    
    mu_prop<-drop(m*exp(X%*%betaprop)/(1+exp(X%*%betaprop)))
    W_prop<-mu_prop*(m-mu_prop)/m
    yhat_prop<-yhat(y,X,betaprop,W_prop,mu_prop)
    
    B.prop<-solve(Binv+t(X)%*%(W_prop*X))
    b.prop<-drop(B.%*%(Binv%*%b+t(X)%*%(W_prop*yhat_prop)))
    
    r<-exp(dmvnorm(t(betaprop),b,B,log=T)-dmvnorm(chain[i-1,],b,B,log=T)+
             sum(dbinom(y,m,mu_prop/m,log =T)-dbinom(y,m,mu_n/m,log =T))+
             dmvnorm(chain[i-1,],b.prop,B.prop,log = T)-(dmvnorm(t(betaprop),b.,B.,log = T)))
    
    if(runif(1)<r){
      accept<-accept+1
      chain[i,]<-betaprop
    } else{
      chain[i,]<-chain[i-1,]
    }
    if(flag){cat(i,accept/N,"\n")}
    Dev[i]<- -2*log_lik(y,m,X,chain[i,])
  }
  return(list(chain=chain,Deviance=Dev,Accepted_samples=accept))
}
binommh.formula <- function(formula, data=list(), weights,...)
{
  mf <- model.frame(formula=formula, data=data)
  X <- model.matrix(attr(mf, "terms"), data=mf)
  Y <- model.response(mf)
  est <- binommh.default(Y,weights, X, ...)
}

poissmh<- function(y,...) UseMethod("poissmh")
poissmh.default<-function(y,X,b=rep(0,dim(X)[2]),B=diag(rep(10000,dim(X)[2])),N=3000,flag=F,...){
  ########utility functions
  yhat<-function(y,X,bta,W_n=1,mu_n){drop(X%*%bta+(y-mu_n)/W_n)}
  btahat<-function(yhat_n,X,W_n=1){drop(solve(t(X)%*%(W_n*X))%*%t(X)%*%(W_n*yhat_n))}
  log_lik<-function(y,X,bta){
    eta_y<-X%*%bta
    sum(dpois(y,drop(exp(eta_y)),log=T)) 
  }
  ####################
  chain<-matrix(rep(NA,dim(X)[2]*N),nrow=N);colnames(chain)<-colnames(X)
  Dev<-rep(NA,N)
  accept <- 0
  ###frcuent cuantities
  Binv<-solve(B)
  ####initials
  ytemp<-y
  ytemp[ytemp==0]<-0.0000001
  W_0<-ytemp
  yhat_0<-log(ytemp)
  chain[1,]<-btahat(yhat_0,X,W_0)
  Dev[1]<- -2*log_lik(y,X,chain[1,])
  ###proposal kernel
  for(i in 2:N){
    mu_n<-drop(exp(X%*%chain[i-1,]))
    W_n<-mu_n
    yhat_n<-yhat(y,X,chain[i-1,],W_n,mu_n)
    
    B.<-solve(Binv+t(X)%*%(W_n*X))
    b.<-drop(B.%*%(Binv%*%b+t(X)%*%(W_n*yhat_n)))
    betaprop<-t(rmvnorm(1,b.,B.))
    
    mu_prop<-drop(exp(X%*%betaprop))
    W_prop<-mu_prop
    yhat_prop<-yhat(y,X,betaprop,W_prop,mu_prop)
    
    B.prop<-solve(Binv+t(X)%*%(W_prop*X))
    b.prop<-drop(B.prop%*%(Binv%*%b+t(X)%*%(W_prop*yhat_prop)))
    r<-exp(dmvnorm(t(betaprop),b,B,log=T)-dmvnorm(chain[i-1,],b,B,log=T)+
             sum(dpois(y,mu_prop,log=T)-dpois(y,mu_n,log=T))+
             dmvnorm(chain[i-1,],b.prop,B.prop,log=T)-dmvnorm(t(betaprop),b.,B.,log=T))
    
    if(runif(1)<r){
      accept<-accept+1
      chain[i,]<-betaprop
    } else{
      chain[i,]<-chain[i-1,]
    }
    if(flag){cat(i,accept/N,"\n")}
    Dev[i]<- -2*log_lik(y,X,chain[i,])
  }
  return(list(chain=chain,Deviance=Dev,Accepted_samples=accept))
}
poissmh.formula <- function(formula, data=list(),...)
{
  mf <- model.frame(formula=formula, data=data)
  X <- model.matrix(attr(mf, "terms"), data=mf)
  Y <- model.response(mf)
  est <- poissmh.default(Y, X, ...)
}

expmh<- function(y,...) UseMethod("expmh")
expmh.default<-function(y,X,b=rep(0,dim(X)[2]),B=diag(rep(100,dim(X)[2])),N=3000,flag=F,...){
  ########utility functions
  yhat<-function(y,X,bta,W_n=1,mu_n){drop(X%*%bta+(y-mu_n)/W_n)}
  btahat<-function(yhat_n,X,W_n=1){drop(solve(t(X)%*%(W_n*X))%*%t(X)%*%(W_n*yhat_n))}
  log_lik<-function(y,X,bta){
    eta_y<-X%*%bta
    sum(dexp(y,1/drop(exp(eta_y)),log=T)) 
  }
  ####################
  chain<-matrix(rep(NA,dim(X)[2]*N),nrow=N);colnames(chain)<-colnames(X)
  Dev<-rep(NA,N)
  accept <- 0
  ###frcuent cuantities
  Binv<-solve(B)
  B.<-solve(Binv+t(X)%*%X)
  ####initials
  ytemp<-y
  ytemp[ytemp==0]<-0.0000001
  yhat_0<-log(ytemp)
  chain[1,]<-btahat(yhat_0,X)
  Dev[1]<- -2*log_lik(y,X,chain[1,])
  ###proposal kernel
  for(i in 2:N){
    mu_n<-drop(exp(X%*%chain[i-1,]))
    W_n<-mu_n
    yhat_n<-yhat(y,X,chain[i-1,],W_n,mu_n)
    
    
    b.<-B.%*%(Binv%*%b+t(X)%*%yhat_n)
    betaprop<-t(rmvnorm(1,b.,B.))
    
    mu_prop<-drop(exp(X%*%betaprop))
    W_prop<-mu_prop
    yhat_prop<-yhat(y,X,betaprop,W_prop,mu_prop)
    
    b.prop<-B.%*%(Binv%*%b+t(X)%*%yhat_prop)
    r<-exp(dmvnorm(t(betaprop),b,B,log=T)-dmvnorm(chain[i-1,],b,B,log=T)+
             sum(dexp(y,1/mu_prop,log =T)-dexp(y,1/mu_n,log =T))+
             dmvnorm(chain[i-1,],b.prop,B.,log = T)-(dmvnorm(t(betaprop),b.,B.,log = T)))
    
    if(runif(1)<r){
      accept<-accept+1
      chain[i,]<-betaprop
    } else{
      chain[i,]<-chain[i-1,]
    }
    Dev[i]<- -2*log_lik(y,X,chain[i,])
    if(flag){cat(i,accept/N,"\n")}
  }
  return(list(chain=chain,Deviance=Dev,Accepted_samples=accept)) 
}
expmh.formula <- function(formula, data=list(),...)
{
  mf <- model.frame(formula=formula, data=data)
  X <- model.matrix(attr(mf, "terms"), data=mf)
  Y <- model.response(mf)
  est <- expmh.default(Y, X, ...)
}

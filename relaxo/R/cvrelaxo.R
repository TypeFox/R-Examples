`cvrelaxo` <-
function(X, Y, K=5, phi= seq(0,1, length=10),  max.steps= min(2*length(Y),2*ncol(X)), fast = TRUE, keep.data=TRUE, warn=TRUE){
Y <- as.numeric(Y)
  if(warn){
    if( abs(mean(Y))> 0.01*sd(Y)) warning("response variable not centered")
    if( any( abs(apply(X,2,mean)) > 0.01* apply(X,2,sd) )) warning("predictor variables not centered")
    if( sd(as.numeric(apply(X,2,sd)))>0.001) warning("predictor variables not scaled")
  }

  n <- length(Y)
  index <- sample( rep(1:K,each=ceiling(n/K)), n, replace=FALSE ) 

  losscv <- rep(0,length=length(phi)*(max.steps-1))
  for (k in 1:K){
    rel <- relaxo(X[index!=k,],Y[index!=k],phi=phi,fast=fast,keep.data=FALSE,warn=FALSE,max.steps=max.steps)
    
    pred <-  X[index==k,]%*%t(rel$beta)
    losscv[1:ncol(pred)] <- losscv[1:ncol(pred)] + apply( sweep(pred,1,Y[index==k])^2,2, mean )/K
    if(length(losscv)>ncol(pred))  losscv[ (ncol(pred)+1):length(losscv)] <-  Inf
  }

  relall <- relaxo(X,Y,phi=phi,fast=fast,keep.data=keep.data,warn=FALSE)
  select <- which.min(losscv)

  relall$beta <- relall$beta[select, ,drop=FALSE]
  relall$lambda <- relall$lambda[select]
  relall$phi <- relall$phi[select]
  
  
  return(relall)
}


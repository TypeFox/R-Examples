CBIV <- function(Tr, Z, X, pZ, method="over", iterations=NULL, standardize = TRUE, twostep = TRUE) {
  probs.min<-10^-6
  
  k<-0
  
  score.only<-bal.only<-FALSE
  if(method=="mle") score.only<-TRUE
  if(method=="exact") bal.only<-TRUE
                       
  X<-cbind(1,X[,apply(X,2,sd)>0])
  names.X<-colnames(X)
  names.X[apply(X,2,sd)==0]<-"(Intercept)"
  
  #######Declare some constants and orthogonalize Xdf.
  X.orig<-X
  format.bal<-F
  x.sd<-apply(as.matrix(X[,-1]),2,sd)
  Dx.inv<-diag(c(1,x.sd))
  diag(Dx.inv)<-1
  x.mean<-apply(as.matrix(X[,-1]),2,mean)
  X[,-1]<-apply(as.matrix(X[,-1]),2,FUN=function(x) (x-mean(x))/sd(x))
  if(k==0) k<-sum(diag(t(X)%*%X%*%ginv(t(X)%*%X)))
  k<-floor(k+.1)
  svd1<-svd(X)
  X<-svd1$u[,1:k]
  XprimeX.inv<-ginv(t(X)%*%X)
  
  if (is.null(iterations)) iterations<-1000
  n<-length(Tr)
  
  gmm.func <- function(beta.curr, invV = NULL)
  {
    beta.curr.c<-beta.curr[1:k]
    beta.curr.a<-beta.curr[k+(1:k)]
    
    baseline.prob<-(1 + exp(X%*%beta.curr.c) + exp(X%*%beta.curr.a))^-1
    
    probs.curr.c<-pmin(pmax(exp(X%*%beta.curr.c)*baseline.prob,probs.min),1-probs.min)
    probs.curr.a<-pmin(pmax(exp(X%*%beta.curr.a)*baseline.prob,probs.min),1-probs.min)
    probs.curr.n<-pmin(pmax(baseline.prob,probs.min),1-probs.min)
    
    sums<-probs.curr.c+probs.curr.a+probs.curr.n
    probs.curr.c<-probs.curr.c/sums
    probs.curr.a<-probs.curr.a/sums
    probs.curr.n<-probs.curr.n/sums
    
    w.curr<-cbind(Z*Tr/(pZ*(probs.curr.c + probs.curr.a)) + (1-Z)*Tr/((1-pZ)*probs.curr.a) - Z*(1-Tr)/(pZ*probs.curr.n) - (1-Z)*(1-Tr)/((1-pZ)*(probs.curr.c + probs.curr.n)),
                  Z*Tr/(pZ*(probs.curr.c + probs.curr.a)) - (1-Z)*Tr/((1-pZ)*probs.curr.a) + Z*(1-Tr)/(pZ*probs.curr.n) - (1-Z)*(1-Tr)/((1-pZ)*(probs.curr.c + probs.curr.n)),
                  Z*Tr/(pZ*(probs.curr.c + probs.curr.a)) - (1-Z)*Tr/((1-pZ)*probs.curr.a) - Z*(1-Tr)/(pZ*probs.curr.n) + (1-Z)*(1-Tr)/((1-pZ)*(probs.curr.c + probs.curr.n)))
    
    w.curr.del<-1/n*t(X)%*%w.curr
    w.curr.del<-as.matrix(w.curr.del)
    w.curr<-as.matrix(w.curr)
    
    gbar<-c(1/n*t(X)%*%((Z*Tr/(probs.curr.c + probs.curr.a) + (1-Z)*(1-Tr)/(1 - probs.curr.a) - 1)*probs.curr.c),
            1/n*t(X)%*%((Z*Tr/(probs.curr.c + probs.curr.a) + (1-Z)*Tr/probs.curr.a - 1)*probs.curr.a),
            w.curr.del)
    
    if (is.null(invV))
    {
      X.1.1<-X*as.vector((pZ/(probs.curr.c + probs.curr.a) + (1 - pZ)/(1 - probs.curr.a) - 1)*probs.curr.c^2)
      X.1.2<-X*as.vector((pZ/(probs.curr.c + probs.curr.a) - 1)*probs.curr.a*probs.curr.c)
      X.1.3<-X*as.vector(probs.curr.c*((probs.curr.c+probs.curr.a)^-1 - (1-probs.curr.a)^-1))
      X.1.4<-X*as.vector(probs.curr.c*((probs.curr.c+probs.curr.a)^-1 - (1-probs.curr.a)^-1))
      X.1.5<-X*as.vector(probs.curr.c*((probs.curr.c+probs.curr.a)^-1 + (1-probs.curr.a)^-1))
      X.2.2<-X*as.vector((pZ/(probs.curr.c + probs.curr.a) + (1 - pZ)/probs.curr.a - 1)*probs.curr.a^2)
      X.2.3<-X*as.vector(probs.curr.a*((probs.curr.c + probs.curr.a)^-1 + probs.curr.a^-1))
      X.2.4<-X*as.vector(probs.curr.a*((probs.curr.c + probs.curr.a)^-1 - probs.curr.a^-1))
      X.2.5<-X*as.vector(probs.curr.a*((probs.curr.c + probs.curr.a)^-1 - probs.curr.a^-1))
      X.3.3<-X*as.vector((pZ*(probs.curr.c + probs.curr.a))^-1 + (pZ*probs.curr.n)^-1 + ((1-pZ)*probs.curr.a)^-1 + ((1-pZ)*(1-probs.curr.a))^-1)
      X.3.4<-X*as.vector((pZ*(probs.curr.c + probs.curr.a))^-1 - (pZ*probs.curr.n)^-1 - ((1-pZ)*probs.curr.a)^-1 + ((1-pZ)*(1-probs.curr.a))^-1)
      X.3.5<-X*as.vector((pZ*(probs.curr.c + probs.curr.a))^-1 + (pZ*probs.curr.n)^-1 - ((1-pZ)*probs.curr.a)^-1 - ((1-pZ)*(1-probs.curr.a))^-1)
      X.4.4<-X*as.vector((pZ*(probs.curr.c + probs.curr.a))^-1 + (pZ*probs.curr.n)^-1 + ((1-pZ)*probs.curr.a)^-1 + ((1-pZ)*(1-probs.curr.a))^-1)
      X.4.5<-X*as.vector((pZ*(probs.curr.c + probs.curr.a))^-1 - (pZ*probs.curr.n)^-1 + ((1-pZ)*probs.curr.a)^-1 - ((1-pZ)*(1-probs.curr.a))^-1)
      X.5.5<-X*as.vector((pZ*(probs.curr.c + probs.curr.a))^-1 + (pZ*probs.curr.n)^-1 + ((1-pZ)*probs.curr.a)^-1 + ((1-pZ)*(1-probs.curr.a))^-1)
      
      V<-1/n*rbind(cbind(t(X.1.1)%*%X.1.1, t(X.1.2)%*%X.1.2, t(X.1.3)%*%X.1.3, t(X.1.4)%*%X.1.4, t(X.1.5)%*%X.1.5),
                   cbind(t(X.1.2)%*%X.1.2, t(X.2.2)%*%X.2.2, t(X.2.3)%*%X.2.3, t(X.2.4)%*%X.2.4, t(X.2.5)%*%X.2.5),
                   cbind(t(X.1.3)%*%X.1.3, t(X.2.3)%*%X.2.3, t(X.3.3)%*%X.3.3, t(X.3.4)%*%X.3.4, t(X.3.5)%*%X.3.5),
                   cbind(t(X.1.4)%*%X.1.4, t(X.2.4)%*%X.2.4, t(X.3.4)%*%X.3.4, t(X.4.4)%*%X.4.4, t(X.4.5)%*%X.4.5),
                   cbind(t(X.1.5)%*%X.1.5, t(X.2.5)%*%X.2.5, t(X.3.5)%*%X.3.4, t(X.4.5)%*%X.4.5, t(X.5.5)%*%X.5.5))
      invV<-ginv(V)
    }
        
    loss1<-as.vector(t(gbar)%*%invV%*%(gbar))      
    out1<-list("loss"=loss1*n, "invV"=invV)
    out1
  }
  
  gmm.loss <- function(beta.curr, invV = NULL) gmm.func(beta.curr, invV)$loss
  
  gmm.gradient <- function(beta.curr, invV)
  {
    beta.curr.c<-beta.curr[1:k]
    beta.curr.a<-beta.curr[k+(1:k)]
    
    baseline.prob<-(1 + exp(X%*%beta.curr.c) + exp(X%*%beta.curr.a))^-1
    
    probs.curr.c<-pmin(pmax(exp(X%*%beta.curr.c)*baseline.prob,probs.min),1-probs.min)
    probs.curr.a<-pmin(pmax(exp(X%*%beta.curr.a)*baseline.prob,probs.min),1-probs.min)
    probs.curr.n<-pmin(pmax(baseline.prob,probs.min),1-probs.min)
    
    sums<-probs.curr.c+probs.curr.a+probs.curr.n
    probs.curr.c<-probs.curr.c/sums
    probs.curr.a<-probs.curr.a/sums
    probs.curr.n<-probs.curr.n/sums
    
    w.curr<-cbind(Z*Tr/(pZ*(probs.curr.c + probs.curr.a)) + (1-Z)*Tr/((1-pZ)*probs.curr.a) - Z*(1-Tr)/(pZ*probs.curr.n) - (1-Z)*(1-Tr)/((1-pZ)*(probs.curr.c + probs.curr.n)),
                  Z*Tr/(pZ*(probs.curr.c + probs.curr.a)) - (1-Z)*Tr/((1-pZ)*probs.curr.a) + Z*(1-Tr)/(pZ*probs.curr.n) - (1-Z)*(1-Tr)/((1-pZ)*(probs.curr.c + probs.curr.n)),
                  Z*Tr/(pZ*(probs.curr.c + probs.curr.a)) - (1-Z)*Tr/((1-pZ)*probs.curr.a) - Z*(1-Tr)/(pZ*probs.curr.n) + (1-Z)*(1-Tr)/((1-pZ)*(probs.curr.c + probs.curr.n)))
    
    w.curr.del<-1/n*t(X)%*%w.curr
    w.curr.del<-as.matrix(w.curr.del)
    w.curr<-as.matrix(w.curr)
    
    gbar<-c(1/n*t(X)%*%((Z*Tr/(probs.curr.c + probs.curr.a) + (1-Z)*(1-Tr)/(1 - probs.curr.a) - 1)*probs.curr.c),
            1/n*t(X)%*%((Z*Tr/(probs.curr.c + probs.curr.a) + (1-Z)*Tr/probs.curr.a - 1)*probs.curr.a),
            w.curr.del)
    
    Ac<- -probs.curr.c*probs.curr.n/(probs.curr.c + probs.curr.a)^2
    Bc<- probs.curr.c/probs.curr.a
    Cc<- -probs.curr.c*probs.curr.a/(1-probs.curr.a)^2
    Dc<- probs.curr.c/probs.curr.n
    Aa<- -probs.curr.a*probs.curr.n/(probs.curr.c + probs.curr.a)^2
    Ba<- -(1-probs.curr.a)/probs.curr.a
    Ca<- probs.curr.a/(1 - probs.curr.a)
    Da<- probs.curr.a/probs.curr.n
    
    dgbar<-rbind(cbind(t(X*as.vector(probs.curr.c*(Z*Tr*Ac + (1-Z)*(1-Tr)*Cc + 
                        (Z*Tr/(probs.curr.c + probs.curr.a) + (1-Z)*(1-Tr)/(1-probs.curr.a) - 1)*(1 - probs.curr.c))))%*%X,
                       t(X*as.vector(probs.curr.a*(Z*Tr*Ac + (1-Z)*Tr*Bc - 
                        (Z*Tr/(probs.curr.c + probs.curr.a) + (1-Z)*Tr/probs.curr.a - 1)*probs.curr.c)))%*%X,
                       t(X*as.vector(Z*Tr/pZ*Ac - Z*(1-Tr)/pZ*Dc + (1-Z)*Tr/(1-pZ)*Bc - 
                                       (1-Z)*(1-Tr)/(1-pZ)*Cc))%*%X,
                       t(X*as.vector(Z*Tr/pZ*Ac + Z*(1-Tr)/pZ*Dc - (1-Z)*Tr/(1-pZ)*Bc - 
                                       (1-Z)*(1-Tr)/(1-pZ)*Cc))%*%X,
                       t(X*as.vector(Z*Tr/pZ*Ac - Z*(1-Tr)/pZ*Dc - (1-Z)*Tr/(1-pZ)*Bc + 
                                       (1-Z)*(1-Tr)/(1-pZ)*Cc))%*%X),
                 cbind(t(X*as.vector(probs.curr.c*(Z*Tr*Aa + (1-Z)*(1-Tr)*Ca - 
                        (Z*Tr/(probs.curr.c + probs.curr.a) + (1-Z)*(1-Tr)/(1-probs.curr.a) - 1)*probs.curr.a)))%*%X,
                       t(X*as.vector(probs.curr.a*(Z*Tr*Aa + (1-Z)*Tr*Ba + 
                        (Z*Tr/(probs.curr.c + probs.curr.a) + (1-Z)*Tr/probs.curr.a - 1)*(1-probs.curr.a))))%*%X,
                       t(X*as.vector(Z*Tr/pZ*Aa - Z*(1-Tr)/pZ*Da + (1-Z)*Tr/(1-pZ)*Ba - 
                                       (1-Z)*(1-Tr)/(1-pZ)*Ca))%*%X,
                       t(X*as.vector(Z*Tr/pZ*Aa + Z*(1-Tr)/pZ*Da - (1-Z)*Tr/(1-pZ)*Ba - 
                                       (1-Z)*(1-Tr)/(1-pZ)*Ca))%*%X,
                       t(X*as.vector(Z*Tr/pZ*Aa - Z*(1-Tr)/pZ*Da - (1-Z)*Tr/(1-pZ)*Ba + 
                                       (1-Z)*(1-Tr)/(1-pZ)*Ca))%*%X))

    out<-2*dgbar%*%invV%*%gbar
    out
  }
  
  mle.loss <- function(beta.curr)
  {
    beta.curr.c<-beta.curr[1:k]
    beta.curr.a<-beta.curr[k+(1:k)]
    
    baseline.prob<-(1 + exp(X%*%beta.curr.c) + exp(X%*%beta.curr.a))^-1
    
    probs.curr.c<-pmin(pmax(exp(X%*%beta.curr.c)*baseline.prob,probs.min),1-probs.min)
    probs.curr.a<-pmin(pmax(exp(X%*%beta.curr.a)*baseline.prob,probs.min),1-probs.min)
    probs.curr.n<-pmin(pmax(1-probs.curr.c-probs.curr.a,probs.min),1-probs.min)
    
    sums<-probs.curr.c+probs.curr.a+probs.curr.n
    probs.curr.c<-probs.curr.c/sums
    probs.curr.a<-probs.curr.a/sums
    probs.curr.n<-probs.curr.n/sums
    
    loss<-1/n*sum(Z*Tr*log(probs.curr.c+probs.curr.a) + Z*(1-Tr)*log(probs.curr.n) + (1-Z)*Tr*log(probs.curr.a) + (1-Z)*(1-Tr)*log(1-probs.curr.a))
    loss
  }
  
  mle.gradient <- function(beta.curr)
  {
    beta.curr.c<-beta.curr[1:k]
    beta.curr.a<-beta.curr[k+(1:k)]
    
    baseline.prob<-(1 + exp(X%*%beta.curr.c) + exp(X%*%beta.curr.a))^-1
    
    probs.curr.c<-pmin(pmax(exp(X%*%beta.curr.c)*baseline.prob,probs.min),1-probs.min)
    probs.curr.a<-pmin(pmax(exp(X%*%beta.curr.a)*baseline.prob,probs.min),1-probs.min)
    probs.curr.n<-pmin(pmax(baseline.prob,probs.min),1-probs.min)
    
    sums<-probs.curr.c+probs.curr.a+probs.curr.n
    probs.curr.c<-probs.curr.c/sums
    probs.curr.a<-probs.curr.a/sums
    probs.curr.n<-probs.curr.n/sums
    
    ds<-1/n*c(t(X)%*%((Z*Tr/(probs.curr.c + probs.curr.a) + (1-Z)*(1-Tr)/(1 - probs.curr.a) - 1)*probs.curr.c),
              t(X)%*%((Z*Tr/(probs.curr.c + probs.curr.a) + (1-Z)*Tr/probs.curr.a - 1)*probs.curr.a))
    ds
  }
  
  bal.loss <- function(beta.curr)
  {
    beta.curr.c<-beta.curr[1:k]
    beta.curr.a<-beta.curr[k+(1:k)]
    
    baseline.prob<-(1 + exp(X%*%beta.curr.c) + exp(X%*%beta.curr.a))^-1
    
    probs.curr.c<-pmin(pmax(exp(X%*%beta.curr.c)*baseline.prob,probs.min),1-probs.min)
    probs.curr.a<-pmin(pmax(exp(X%*%beta.curr.a)*baseline.prob,probs.min),1-probs.min)
    probs.curr.n<-pmin(pmax(baseline.prob,probs.min),1-probs.min)
    
    sums<-probs.curr.c+probs.curr.a+probs.curr.n
    probs.curr.c<-probs.curr.c/sums
    probs.curr.a<-probs.curr.a/sums
    probs.curr.n<-probs.curr.n/sums
    
    w.curr<-cbind(Z*Tr/(pZ*(probs.curr.c + probs.curr.a)) + (1-Z)*Tr/((1-pZ)*probs.curr.a) - Z*(1-Tr)/(pZ*probs.curr.n) - (1-Z)*(1-Tr)/((1-pZ)*(probs.curr.c + probs.curr.n)),
                  Z*Tr/(pZ*(probs.curr.c + probs.curr.a)) - (1-Z)*Tr/((1-pZ)*probs.curr.a) + Z*(1-Tr)/(pZ*probs.curr.n) - (1-Z)*(1-Tr)/((1-pZ)*(probs.curr.c + probs.curr.n)),
                  Z*Tr/(pZ*(probs.curr.c + probs.curr.a)) - (1-Z)*Tr/((1-pZ)*probs.curr.a) - Z*(1-Tr)/(pZ*probs.curr.n) + (1-Z)*(1-Tr)/((1-pZ)*(probs.curr.c + probs.curr.n)))
    
    loss1<-sum(diag(abs(t(w.curr)%*%X%*%XprimeX.inv%*%t(X)%*%(w.curr))))
    loss1
  }
  
  bal.gradient <- function(beta.curr)
  {
    beta.curr.c<-beta.curr[1:k]
    beta.curr.a<-beta.curr[k+(1:k)]
    
    baseline.prob<-(1 + exp(X%*%beta.curr.c) + exp(X%*%beta.curr.a))^-1
    
    probs.curr.c<-pmin(pmax(exp(X%*%beta.curr.c)*baseline.prob,probs.min),1-probs.min)
    probs.curr.a<-pmin(pmax(exp(X%*%beta.curr.a)*baseline.prob,probs.min),1-probs.min)
    probs.curr.n<-pmin(pmax(baseline.prob,probs.min),1-probs.min)
    
    sums<-probs.curr.c+probs.curr.a+probs.curr.n
    probs.curr.c<-probs.curr.c/sums
    probs.curr.a<-probs.curr.a/sums
    probs.curr.n<-probs.curr.n/sums
    
    Ac<- -probs.curr.c*probs.curr.n/(probs.curr.c + probs.curr.a)^2
    Bc<- probs.curr.c/probs.curr.a
    Cc<- -probs.curr.c*probs.curr.a/(1-probs.curr.a)^2
    Dc<- probs.curr.c/probs.curr.n
    Aa<- -probs.curr.a*probs.curr.n/(probs.curr.c + probs.curr.a)^2
    Ba<- -(1-probs.curr.a)/probs.curr.a
    Ca<- probs.curr.a/(1 - probs.curr.a)
    Da<- probs.curr.a/probs.curr.n
    
    w.curr<-cbind(Z*Tr/(pZ*(probs.curr.c + probs.curr.a)) + (1-Z)*Tr/((1-pZ)*probs.curr.a) - Z*(1-Tr)/(pZ*probs.curr.n) - (1-Z)*(1-Tr)/((1-pZ)*(probs.curr.c + probs.curr.n)),
                  Z*Tr/(pZ*(probs.curr.c + probs.curr.a)) - (1-Z)*Tr/((1-pZ)*probs.curr.a) + Z*(1-Tr)/(pZ*probs.curr.n) - (1-Z)*(1-Tr)/((1-pZ)*(probs.curr.c + probs.curr.n)),
                  Z*Tr/(pZ*(probs.curr.c + probs.curr.a)) - (1-Z)*Tr/((1-pZ)*probs.curr.a) - Z*(1-Tr)/(pZ*probs.curr.n) + (1-Z)*(1-Tr)/((1-pZ)*(probs.curr.c + probs.curr.n)))
    
    dw.beta.c<-cbind(t(X*as.vector(Z*Tr/pZ*Ac - Z*(1-Tr)/pZ*Dc + (1-Z)*Tr/(1-pZ)*Bc - (1-Z)*(1-Tr)/(1-pZ)*Cc)),
                     t(X*as.vector(Z*Tr/pZ*Ac + Z*(1-Tr)/pZ*Dc - (1-Z)*Tr/(1-pZ)*Bc - (1-Z)*(1-Tr)/(1-pZ)*Cc)),
                     t(X*as.vector(Z*Tr/pZ*Ac - Z*(1-Tr)/pZ*Dc - (1-Z)*Tr/(1-pZ)*Bc + (1-Z)*(1-Tr)/(1-pZ)*Cc)))

    dw.beta.a<-cbind(t(X*as.vector(Z*Tr/pZ*Aa - Z*(1-Tr)/pZ*Da + (1-Z)*Tr/(1-pZ)*Ba - (1-Z)*(1-Tr)/(1-pZ)*Ca)),
                     t(X*as.vector(Z*Tr/pZ*Aa + Z*(1-Tr)/pZ*Da - (1-Z)*Tr/(1-pZ)*Ba - (1-Z)*(1-Tr)/(1-pZ)*Ca)),
                     t(X*as.vector(Z*Tr/pZ*Aa - Z*(1-Tr)/pZ*Da - (1-Z)*Tr/(1-pZ)*Ba + (1-Z)*(1-Tr)/(1-pZ)*Ca)))

    out.1<-2*dw.beta.c[,1:n]%*%X%*%XprimeX.inv%*%t(X)%*%(w.curr[,1]) + 
           2*dw.beta.c[,n+(1:n)]%*%X%*%XprimeX.inv%*%t(X)%*%(w.curr[,2]) +
           2*dw.beta.c[,2*n+(1:n)]%*%X%*%XprimeX.inv%*%t(X)%*%(w.curr[,3])
    out.2<-2*dw.beta.a[,1:n]%*%X%*%XprimeX.inv%*%t(X)%*%(w.curr[,1]) + 
           2*dw.beta.a[,n+(1:n)]%*%X%*%XprimeX.inv%*%t(X)%*%(w.curr[,2]) +
           2*dw.beta.a[,2*n+(1:n)]%*%X%*%XprimeX.inv%*%t(X)%*%(w.curr[,3])
  
    out<-c(out.1, out.2)
    out
  }
  
  beta.init<-rep(0,2*k)
  mle.opt<-optim(beta.init, mle.loss, control=list("maxit"=iterations), method = "BFGS", gr = mle.gradient)
  beta.mle<-mle.opt$par
    
  if (score.only)   gmm.opt<-mle.opt
  else {
    bal.opt<-optim(beta.mle, bal.loss, control=list("maxit"=iterations), method = "BFGS", gr = bal.gradient)
    beta.bal<-bal.opt$par
    
    #print(bal.gradient(beta.bal))
    #print(grad(bal.loss,beta.bal))
    
    this.invV<-gmm.func(beta.mle)$invV
    
    if (bal.only) gmm.opt<-bal.opt
    else {
      gmm.mle.opt<-optim(beta.mle, gmm.loss, control=list("maxit"=iterations), method = "BFGS", invV = this.invV, gr = gmm.gradient)
      gmm.bal.opt<-optim(beta.bal, gmm.loss, control=list("maxit"=iterations), method = "BFGS", invV = this.invV, gr = gmm.gradient)      
      if (gmm.mle.opt$value > gmm.bal.opt$value)
      {
        gmm.opt<-gmm.bal.opt
      }
      else
      {
        gmm.opt<-gmm.mle.opt
      }
      
      #print(gmm.gradient(gmm.opt$par, invV=this.invV))
      #print(grad(function(x) gmm.loss(x, invV=this.invV),gmm.opt$par))
    }
  }

  #print(gmm.opt$par)
  J.opt<-matrix(gmm.opt$val)
  beta.opt<-matrix(gmm.opt$par,nrow=k)
  class(beta.opt)<-"coef"
  
  baseline<-(1+exp(X%*%beta.opt[,1])+exp(X%*%beta.opt[,2]))^-1
  fitted.values<-cbind(pmin(pmax(exp(X%*%beta.opt[,1]),probs.min),1-probs.min),pmin(pmax(exp(X%*%beta.opt[,2]),probs.min),1-probs.min))
  fitted.values<-cbind(fitted.values, apply(fitted.values, 1, function(x) min(max(1 - sum(x),probs.min),1-probs.min)))
  sums<-apply(fitted.values,1,sum)
  fitted.values[,1]<-fitted.values[,1]/sums
  fitted.values[,2]<-fitted.values[,2]/sums
  fitted.values[,3]<-fitted.values[,3]/sums
  colnames(fitted.values)<-c("Compliers","Always","Never")  
  
  d.inv<- svd1$d
  d.inv[d.inv> 1e-5]<-1/d.inv[d.inv> 1e-5]
  d.inv[d.inv<= 1e-5]<-0
  beta.opt<-svd1$v%*%diag(d.inv)%*%beta.opt
  beta.opt[-1,]<-beta.opt[-1,]/x.sd
  
  deviance<- -2*sum(Z*Tr*log(fitted.values[,1]+fitted.values[,2]) + Z*(1-Tr)*log(fitted.values[,3]) + (1-Z)*Tr*log(fitted.values[,2]) + (1-Z)*(1-Tr)*log(1-fitted.values[,2]))
  
  if (k > 2)
  {
    beta.opt[1,]<-beta.opt[1,]-matrix(x.mean%*%beta.opt[-1,])
  }
  else
  {
    beta.opt[1,]<-beta.opt[1,]-x.mean*beta.opt[-1,]
  }

  output<-list("coefficients"=beta.opt,"fitted.values"=X%*%beta.opt,"rank"=k,
               "deviance"=deviance,"converged"=gmm.opt$conv,"J"=J.opt,"df"=k,
               "bal"=bal.loss(beta.opt))
  class(output)<-"CBIV"
  output
}
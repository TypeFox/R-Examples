CBPS.3Treat<-function(treat, X, X.bal, method, k, XprimeX.inv, bal.only, iterations, standardize, twostep, ...)
{
  probs.min<-1e-6
  
  no.treats<-length(levels(as.factor(treat)))
  treat.names<-levels(as.factor(treat))
  T1<-as.numeric(treat==treat.names[1])
  T2<-as.numeric(treat==treat.names[2])
  T3<-as.numeric(treat==treat.names[3])
  
  ##The gmm objective function--given a guess of beta, constructs the GMM J statistic.
  gmm.func<-function(beta.curr,X.gmm=X,invV=NULL){
    ##Designate a few objects in the function.
    beta.curr<-matrix(beta.curr,k,no.treats-1)
    X<-as.matrix(X.gmm)
    
    ##Designate sample size, number of treated and control observations,
    ##theta.curr, which are used to generate probabilities.
    ##Trim probabilities, and generate weights.
    n<-c(sum(T1), sum(T2), sum(T3))
    theta.curr<-X%*%beta.curr
    baseline.prob<-apply(theta.curr,1,function(x) (1+sum(exp(x)))^-1)
    probs.curr<-cbind(baseline.prob, exp(theta.curr[,1])*baseline.prob, exp(theta.curr[,2])*baseline.prob)
    probs.curr[,1]<-pmax(probs.min,probs.curr[,1])
    probs.curr[,2]<-pmax(probs.min,probs.curr[,2])
    probs.curr[,3]<-pmax(probs.min,probs.curr[,3])
    norms<-apply(probs.curr,1,sum)
    probs.curr<-probs.curr/norms
    
    w.curr<-cbind(2*T1/probs.curr[,1] - T2/probs.curr[,2] - T3/probs.curr[,3],
                  T2/probs.curr[,2] - T3/probs.curr[,3])
    
    
    ##Generate the vector of mean imbalance by weights.
    w.curr.del<-1/sum(n)*t(X)%*%w.curr
    w.curr.del<-as.matrix(w.curr.del)
    w.curr<-as.matrix(w.curr)
    
    ##Generate g-bar, as in the paper.
    gbar<-c(1/sum(n)*t(X)%*%(T2-probs.curr[,2]), 
            1/sum(n)*t(X)%*%(T3-probs.curr[,3]),
            w.curr.del)
    
    
    if(is.null(invV))
    {
      ##Generate the covariance matrix used in the GMM estimate.
      ##Was for the initial version that calculates the analytic variances.    
      X.1.1<-X*(probs.curr[,2]*(1-probs.curr[,2]))
      X.1.2<-X*(-probs.curr[,2]*probs.curr[,3])
      X.1.3<-X*-1
      X.1.4<-X*1
      X.2.2<-X*(probs.curr[,3]*(1-probs.curr[,3]))
      X.2.3<-X*-1
      X.2.4<-X*-1
      X.3.3<-X*(4*probs.curr[,1]^-1+probs.curr[,2]^-1+probs.curr[,3]^-1)
      X.3.4<-X*(-probs.curr[,2]^-1+probs.curr[,3]^-1)
      X.4.4<-X*(probs.curr[,2]^-1+probs.curr[,3]^-1)
      
      V<-1/sum(n)*rbind(cbind(t(X.1.1)%*%X,t(X.1.2)%*%X,t(X.1.3)%*%X,t(X.1.4)%*%X),
                        cbind(t(X.1.2)%*%X,t(X.2.2)%*%X,t(X.2.3)%*%X,t(X.2.4)%*%X),
                        cbind(t(X.1.3)%*%X,t(X.2.3)%*%X,t(X.3.3)%*%X,t(X.3.4)%*%X),
                        cbind(t(X.1.4)%*%X,t(X.2.4)%*%X,t(X.3.4)%*%X,t(X.4.4)%*%X))
      
      invV<-ginv(V)
    }
    
    ##Calculate the GMM loss.
    loss1<-as.vector(t(gbar)%*%invV%*%(gbar))      
    out1<-list("loss"=max(loss1*sum(n),loss1*sum(n)), "invV"=invV)
    out1
  }
  gmm.loss<-function(x,...) gmm.func(x,...)$loss
  
  ##Loss function for balance constraints, returns the squared imbalance along each dimension.
  bal.loss<-function(beta.curr){
    beta.curr<-matrix(beta.curr,k,no.treats-1)
    theta.curr<-X%*%beta.curr
    baseline.prob<-apply(theta.curr,1,function(x) (1+sum(exp(x)))^-1)
    probs.curr<-cbind(baseline.prob, exp(theta.curr[,1])*baseline.prob, exp(theta.curr[,2])*baseline.prob)
    probs.curr[,1]<-pmax(probs.min,probs.curr[,1])
    probs.curr[,2]<-pmax(probs.min,probs.curr[,2])
    probs.curr[,3]<-pmax(probs.min,probs.curr[,3])
    norms<-apply(probs.curr,1,sum)
    probs.curr<-probs.curr/norms
    
    w.curr<-cbind(2*T1/probs.curr[,1] - T2/probs.curr[,2] - T3/probs.curr[,3],
                  T2/probs.curr[,2] - T3/probs.curr[,3])
    ##Generate mean imbalance.
    loss1<-sum(diag(abs(t(w.curr)%*%X%*%XprimeX.inv%*%t(X)%*%(w.curr))))
    loss1
  }
  
  gmm.gradient<-function(beta.curr, X.gmm=X, invV)
  {
    ##Designate a few objects in the function.
    beta.curr<-matrix(beta.curr,k,no.treats-1)
    X<-as.matrix(X.gmm)
    
    ##Designate sample size, number of treated and control observations,
    ##theta.curr, which are used to generate probabilities.
    ##Trim probabilities, and generate weights.
    n<-c(sum(T1), sum(T2), sum(T3))
    theta.curr<-X%*%beta.curr
    baseline.prob<-apply(theta.curr,1,function(x) (1+sum(exp(x)))^-1)
    probs.curr<-cbind(baseline.prob, exp(theta.curr[,1])*baseline.prob, exp(theta.curr[,2])*baseline.prob)
    probs.curr[,1]<-pmax(probs.min,probs.curr[,1])
    probs.curr[,2]<-pmax(probs.min,probs.curr[,2])
    probs.curr[,3]<-pmax(probs.min,probs.curr[,3])
    norms<-apply(probs.curr,1,sum)
    probs.curr<-probs.curr/norms
    
    w.curr<-cbind(2*T1/probs.curr[,1] - T2/probs.curr[,2] - T3/probs.curr[,3],
                  T2/probs.curr[,2] - T3/probs.curr[,3])
    
    
    ##Generate the vector of mean imbalance by weights.
    w.curr.del<-1/sum(n)*t(X)%*%w.curr
    w.curr.del<-as.matrix(w.curr.del)
    w.curr<-as.matrix(w.curr)
    
    ##Generate g-bar, as in the paper.
    gbar<-c(1/sum(n)*t(X)%*%(T2-probs.curr[,2]), 
            1/sum(n)*t(X)%*%(T3-probs.curr[,3]),
            w.curr.del)
    
    dgbar<-rbind(cbind(1/sum(n)*t(apply(t(X),1,function(x) -x*probs.curr[,2]*(1-probs.curr[,2])))%*%X,
                       1/sum(n)*t(apply(t(X),1,function(x) x*probs.curr[,2]*probs.curr[,3]))%*%X,
                       1/sum(n)*t(apply(t(X),1,function(x) x*(2*T1*probs.curr[,2]/probs.curr[,1] + T2*(1-probs.curr[,2])/probs.curr[,2] - T3*probs.curr[,2]/probs.curr[,3])))%*%X,
                       1/sum(n)*t(apply(t(X),1,function(x) x*(-T2*(1-probs.curr[,2])/probs.curr[,2] - T3*probs.curr[,2]/probs.curr[,3])))%*%X),
                 cbind(1/sum(n)*t(apply(t(X),1,function(x) x*probs.curr[,2]*probs.curr[,3]))%*%X,
                       1/sum(n)*t(apply(t(X),1,function(x) -x*probs.curr[,3]*(1-probs.curr[,3])))%*%X,
                       1/sum(n)*t(apply(t(X),1,function(x) x*(2*T1*probs.curr[,3]/probs.curr[,1] - T2*probs.curr[,3]/probs.curr[,2] + T3*(1-probs.curr[,3])/probs.curr[,3])))%*%X,
                       1/sum(n)*t(apply(t(X),1,function(x) x*(T2*probs.curr[,3]/probs.curr[,2] + T3*(1-probs.curr[,3])/probs.curr[,3])))%*%X))
                 
    out<-2*sum(n)*dgbar%*%invV%*%gbar
    out
  }
  
  bal.gradient<-function(beta.curr)
  {
    beta.curr<-matrix(beta.curr,k,no.treats-1)
    theta.curr<-X%*%beta.curr
    baseline.prob<-apply(theta.curr,1,function(x) (1+sum(exp(x)))^-1)
    probs.curr<-cbind(baseline.prob, exp(theta.curr[,1])*baseline.prob, exp(theta.curr[,2])*baseline.prob)
    probs.curr[,1]<-pmax(probs.min,probs.curr[,1])
    probs.curr[,2]<-pmax(probs.min,probs.curr[,2])
    probs.curr[,3]<-pmax(probs.min,probs.curr[,3])
    norms<-apply(probs.curr,1,sum)
    probs.curr<-probs.curr/norms
    
    w.curr<-cbind(2*T1/probs.curr[,1] - T2/probs.curr[,2] - T3/probs.curr[,3],
                  T2/probs.curr[,2] - T3/probs.curr[,3])
    
    dw.beta1<-cbind(t(apply(t(X),1,function(x) x*(2*T1*probs.curr[,2]/probs.curr[,1] + T2*(1-probs.curr[,2])/probs.curr[,2] - T3*probs.curr[,2]/probs.curr[,3]))),
                    t(apply(t(X),1,function(x) x*(-T2*(1-probs.curr[,2])/probs.curr[,2] - T3*probs.curr[,2]/probs.curr[,3]))))
    dw.beta2<-cbind(t(apply(t(X),1,function(x) x*(2*T1*probs.curr[,3]/probs.curr[,1] - T2*probs.curr[,3]/probs.curr[,2] + T3*(1-probs.curr[,3])/probs.curr[,3]))),
                    t(apply(t(X),1,function(x) x*(T2*probs.curr[,3]/probs.curr[,2] + T3*(1-probs.curr[,3])/probs.curr[,3]))))
    ##Generate mean imbalance.
    loss1<-diag(t(w.curr)%*%X%*%XprimeX.inv%*%t(X)%*%(w.curr))
    out.1<-sapply(2*dw.beta1[,1:n]%*%X%*%XprimeX.inv%*%t(X)%*%(w.curr[,1]), 
                  function(x) ifelse((x > 0 & loss1[1] > 0) | (x < 0 & loss1[1] < 0),abs(x),-abs(x))) +
      sapply(2*dw.beta1[,(n+1):(2*n)]%*%X%*%XprimeX.inv%*%t(X)%*%(w.curr[,2]),
             function(x) ifelse((x > 0 & loss1[2] > 0) | (x < 0 & loss1[2] < 0),abs(x),-abs(x)))
    out.2<-sapply(2*dw.beta2[,1:n]%*%X%*%XprimeX.inv%*%t(X)%*%(w.curr[,1]), 
                  function(x) ifelse((x > 0 & loss1[1] > 0) | (x < 0 & loss1[1] < 0),abs(x),-abs(x))) +
      sapply(2*dw.beta2[,(n+1):(2*n)]%*%X%*%XprimeX.inv%*%t(X)%*%(w.curr[,2]),
             function(x) ifelse((x > 0 & loss1[2] > 0) | (x < 0 & loss1[2] < 0),abs(x),-abs(x)))
    out<-c(out.1, out.2)
    out
  }
  
  n<-length(treat)
  n.t<-c(sum(T1),sum(T2),sum(T3))
  x.orig<-x<-cbind(as.matrix(X))
  
  ##Run multionmial logit
  dat.dummy<-data.frame(treat=treat,X)
  #Need to generalize for different dimensioned X's
  xnam<- paste0("X", 1:dim(X)[2])
  fmla <- as.formula(paste("as.factor(treat) ~ -1 + ", paste(xnam, collapse= "+")))
  mnl1<-multinom(fmla, data=dat.dummy, trace=FALSE)
  mcoef<-matrix(coef(mnl1),k,no.treats-1,byrow=TRUE)
  mcoef[is.na(mcoef[,1]),1]<-0
  mcoef[is.na(mcoef[,2]),2]<-0
  probs.mnl<-cbind(1/(1+exp(X%*%mcoef[,1])+exp(X%*%mcoef[,2])),
                   exp(X%*%mcoef[,1])/(1+exp(X%*%mcoef[,1])+exp(X%*%mcoef[,2])),
                   exp(X%*%mcoef[,2])/(1+exp(X%*%mcoef[,1])+exp(X%*%mcoef[,2])))
  colnames(probs.mnl)<-c("p1","p2","p3")
  probs.mnl[,1]<-pmax(probs.min,probs.mnl[,1])
  probs.mnl[,2]<-pmax(probs.min,probs.mnl[,2])
  probs.mnl[,3]<-pmax(probs.min,probs.mnl[,3])
  norms<-apply(probs.mnl,1,sum)
  probs.mnl<-probs.mnl/norms
  
  mnl1$fit<-matrix(probs.mnl,nrow=n,ncol=no.treats)
  beta.curr<-mcoef
  beta.curr[is.na(beta.curr)]<-0
  
  alpha.func<-function(alpha) gmm.loss(beta.curr*alpha)
  beta.curr<-beta.curr*optimize(alpha.func,interval=c(.8,1.1))$min
  
  if(twostep)
  {
    this.invV<-gmm.func(beta.curr)$invV
  }
  
  
  ##Generate estimates for balance and CBPSE
  gmm.init<-beta.curr
  
  if (twostep)
  {
    opt.bal<-optim(gmm.init, bal.loss, control=list("maxit"=iterations), method="BFGS", gr = bal.gradient, hessian=TRUE)
  }
  else
  {
    opt.bal<-tryCatch({
      optim(gmm.init, bal.loss, control=list("maxit"=iterations), method="BFGS", hessian=TRUE)
    },
    error = function(err)
    {
      return(optim(gmm.init, bal.loss, control=list("maxit"=iterations), method="Nelder-Mead", hessian=TRUE))
    })	  
  }
  beta.bal<-opt.bal$par
  if(bal.only) opt1<-opt.bal 
  
  if(!bal.only)
  {
    if (twostep)
    {
      gmm.glm.init<-optim(gmm.init, gmm.loss, control=list("maxit"=iterations), method="BFGS", hessian=TRUE, gr = gmm.gradient, invV = this.invV)
      gmm.bal.init<-optim(beta.bal, gmm.loss, control=list("maxit"=iterations), method="BFGS", hessian=TRUE, gr = gmm.gradient, invV = this.invV)
    }
    else
    {
      gmm.glm.init<-tryCatch({
        optim(gmm.init,gmm.loss, control=list("maxit"=iterations), method="BFGS", hessian=TRUE)
      },
      error = function(err)
      {
        return(optim(gmm.init,gmm.loss, control=list("maxit"=iterations), method="Nelder-Mead", hessian=TRUE))
      })
      gmm.bal.init<-tryCatch({
        optim(beta.bal, gmm.loss, control=list("maxit"=iterations), method="BFGS", hessian=TRUE)
      },
      error = function(err)
      {
        return(optim(beta.bal, gmm.loss, control=list("maxit"=iterations), method="Nelder-Mead", hessian=TRUE))
      })
    }
    
    if(gmm.glm.init$val<gmm.bal.init$val) opt1<-gmm.glm.init else opt1<-gmm.bal.init
  }
  
  ##Generate probabilities
  beta.opt<-matrix(opt1$par,nrow=k,ncol=no.treats-1)
  theta.opt<-X%*%beta.opt
  baseline.prob<-apply(theta.opt,1,function(x) (1+sum(exp(x)))^-1)
  probs.opt<-cbind(baseline.prob, exp(theta.opt[,1])*baseline.prob, exp(theta.opt[,2])*baseline.prob)
  probs.opt[,1]<-pmax(probs.min,probs.opt[,1])
  probs.opt[,2]<-pmax(probs.min,probs.opt[,2])
  probs.opt[,3]<-pmax(probs.min,probs.opt[,3])
  norms<-apply(probs.opt,1,sum)
  probs.opt<-probs.opt/norms
  
  
  J.opt<-ifelse(twostep, gmm.func(beta.opt, invV = this.invV)$loss, gmm.func(beta.opt)$loss)
  
  if ((J.opt > gmm.loss(mcoef)) & (bal.loss(beta.opt) > bal.loss(mcoef)))
  {	  
    beta.opt<-mcoef
    probs.opt<-probs.mnl
    J.opt <- gmm.loss(mcoef)
    warning("Optimization failed.  Results returned are for MLE.")
  }
  
  ##Generate weights
  w.opt<-rbind(2*T1/probs.opt[,1] - T2/probs.opt[,2] - T3/probs.opt[,3],
               T2/probs.opt[,2] - T3/probs.opt[,3])
  
  residuals<-cbind(T1-probs.opt[,1],T2-probs.opt[,2],T3-probs.opt[,3])
  deviance <- -2*c(sum(T1*log(probs.opt[,1])+T2*log(probs.opt[,2])+T3*log(probs.opt[,3])))
  nulldeviance <- -2*c(sum(T1*log(mean(T1))+T2*log(mean(T2))+T3*log(mean(T3))))
  
  norm1<-norm2<-norm3<-1
  if (standardize)
  {
    norm1<-sum(T1/probs.opt[,1])
    norm2<-sum(T2/probs.opt[,2])
    norm3<-sum(T3/probs.opt[,3])
  }
  w<-T1/probs.opt[,1]/norm1 + T2/probs.opt[,2]/norm2 + T3/probs.opt[,3]/norm3
  
  if (twostep)
  {
    W<-this.invV
  }
  else
  {
    W<-gmm.func(beta.opt)$invV
  }
  XG.1.1<-t(-X*probs.opt[,2]*(1-probs.opt[,2]))%*%X
  XG.1.2<-t(X*probs.opt[,2]*probs.opt[,3])%*%X
  XG.1.3<-t(X*(2*T1*probs.opt[,2]/probs.opt[,1] + T2*(1-probs.opt[,2])/probs.opt[,2] - T3*probs.opt[,2]/probs.opt[,3]))%*%X
  XG.1.4<-t(X*(T2*(1-probs.opt[,2])/probs.opt[,2] - T3*probs.opt[,2]/probs.opt[,3]))%*%X
  XG.2.1<-t(X*probs.opt[,2]*probs.opt[,3])%*%X
  XG.2.2<-t(-X*probs.opt[,3]*(1-probs.opt[,3]))%*%X
  XG.2.3<-t(X*(2*T1*probs.opt[,3]/probs.opt[,1] - T2*probs.opt[,3]/probs.opt[,2] + T3*(1-probs.opt[,3])/probs.opt[,3]))%*%X
  XG.2.4<-t(X*(T2*probs.opt[,3]/probs.opt[,2] + T3*(1-probs.opt[,3])/probs.opt[,3]))%*%X
  G<-1/n*rbind(cbind(XG.1.1,XG.1.2,XG.1.3,XG.1.4),cbind(XG.2.1,XG.2.2,XG.2.3,XG.2.4))
  
  XW.1<-X*(T2-probs.opt[,2])
  XW.2<-X*(T3-probs.opt[,3])
  XW.3<-X*(2*T1/probs.opt[,1] - T2/probs.opt[,2] - T3/probs.opt[,3])
  XW.4<-X*(T2/probs.opt[,2] - T3/probs.opt[,3])
  W1<-rbind(t(XW.1),t(XW.2),t(XW.3),t(XW.4))
  Omega<-1/n*(W1%*%t(W1))
  vcov<-ginv(G%*%W%*%t(G))%*%G%*%W%*%Omega%*%W%*%t(G)%*%ginv(G%*%W%*%t(G))
  
  colnames(probs.opt)<-treat.names
  
  class(beta.opt) <- "coef"
  
  output<-list("coefficients"=beta.opt,"fitted.values"=probs.opt,"deviance"=deviance,"weights"=w, 
               "y"=treat,"x"=X,"converged"=opt1$conv,"J"=J.opt,"var"=vcov, 
               "mle.J"=ifelse(twostep, gmm.func(mcoef, invV = this.invV)$loss, gmm.loss(mcoef)))
  
  class(output)<- c("CBPS")
  
  output
}

CBPS.4Treat<-function(treat, X, X.bal, method, k, XprimeX.inv, bal.only, iterations, standardize, twostep, ...)
{
  probs.min<-1e-6
  
  no.treats<-length(levels(as.factor(treat)))
  treat.names<-levels(as.factor(treat))
  T1<-as.numeric(treat==treat.names[1])
  T2<-as.numeric(treat==treat.names[2])
  T3<-as.numeric(treat==treat.names[3])
  T4<-as.numeric(treat==treat.names[4])
  
  ##The gmm objective function--given a guess of beta, constructs the GMM J statistic.
  gmm.func<-function(beta.curr,X.gmm=X,invV=NULL){
    ##Designate a few objects in the function.
    beta.curr<-matrix(beta.curr,k,no.treats-1)
    X<-as.matrix(X.gmm)
    
    ##Designate sample size, number of treated and control observations,
    ##theta.curr, which are used to generate probabilities.
    ##Trim probabilities, and generate weights.
    n<-c(sum(T1), sum(T2), sum(T3), sum(T4))
    theta.curr<-X%*%beta.curr
    baseline.prob<-apply(theta.curr,1,function(x) (1+sum(exp(x)))^-1)
    probs.curr<-cbind(baseline.prob, exp(theta.curr[,1])*baseline.prob, exp(theta.curr[,2])*baseline.prob, exp(theta.curr[,3])*baseline.prob)
    probs.curr[,1]<-pmax(probs.min,probs.curr[,1])
    probs.curr[,2]<-pmax(probs.min,probs.curr[,2])
    probs.curr[,3]<-pmax(probs.min,probs.curr[,3])
    probs.curr[,4]<-pmax(probs.min,probs.curr[,4])
    norms<-apply(probs.curr,1,sum)
    probs.curr<-probs.curr/norms
    
    w.curr<-cbind(T1/probs.curr[,1] + T2/probs.curr[,2] - T3/probs.curr[,3] - T4/probs.curr[,4],
                  T1/probs.curr[,1] - T2/probs.curr[,2] - T3/probs.curr[,3] + T4/probs.curr[,4],
                  -T1/probs.curr[,1] + T2/probs.curr[,2] - T3/probs.curr[,3] + T4/probs.curr[,4])    
    
    ##Generate the vector of mean imbalance by weights.
    w.curr.del<-1/sum(n)*t(X)%*%w.curr
    w.curr.del<-as.matrix(w.curr.del)
    w.curr<-as.matrix(w.curr)
    
    ##Generate g-bar, as in the paper.
    gbar<-c(1/sum(n)*t(X)%*%(T2-probs.curr[,2]), 
            1/sum(n)*t(X)%*%(T3-probs.curr[,3]),
            1/sum(n)*t(X)%*%(T4-probs.curr[,4]),
            w.curr.del)
    
    
    if(is.null(invV))
    {
      ##Generate the covariance matrix used in the GMM estimate.
      ##Was for the initial version that calculates the analytic variances.		
      X.1.1<-X*(probs.curr[,2]*(1-probs.curr[,2]))
      X.1.2<-X*(-probs.curr[,2]*probs.curr[,3])
      X.1.3<-X*(-probs.curr[,2]*probs.curr[,4])
      X.1.4<-X
      X.1.5<-X*(-1)
      X.1.6<-X
      X.2.2<-X*(probs.curr[,3]*(1-probs.curr[,3]))
      X.2.3<-X*(-probs.curr[,3]*probs.curr[,4])
      X.2.4<-X*(-1)
      X.2.5<-X*(-1)
      X.2.6<-X*(-1)
      X.3.3<-X*(probs.curr[,4]*(1-probs.curr[,4]))
      X.3.4<-X*(-1)
      X.3.5<-X
      X.3.6<-X
      X.4.4<-X*(probs.curr[,1]^-1+probs.curr[,2]^-1+probs.curr[,3]^-1+probs.curr[,4]^-1)
      X.4.5<-X*(probs.curr[,1]^-1-probs.curr[,2]^-1+probs.curr[,3]^-1-probs.curr[,4]^-1)
      X.4.6<-X*(-probs.curr[,1]^-1+probs.curr[,2]^-1+probs.curr[,3]^-1-probs.curr[,4]^-1)
      X.5.5<-X.4.4
      X.5.6<-X*(-probs.curr[,1]^-1-probs.curr[,2]^-1+probs.curr[,3]^-1+probs.curr[,4]^-1)
      X.6.6<-X.4.4
      
      V<-1/sum(n)*rbind(cbind(t(X.1.1)%*%X,t(X.1.2)%*%X,t(X.1.3)%*%X,t(X.1.4)%*%X,t(X.1.5)%*%X,t(X.1.6)%*%X),
                        cbind(t(X.1.2)%*%X,t(X.2.2)%*%X,t(X.2.3)%*%X,t(X.2.4)%*%X,t(X.2.5)%*%X,t(X.2.6)%*%X),
                        cbind(t(X.1.3)%*%X,t(X.2.3)%*%X,t(X.3.3)%*%X,t(X.3.4)%*%X,t(X.3.5)%*%X,t(X.3.6)%*%X),
                        cbind(t(X.1.4)%*%X,t(X.2.4)%*%X,t(X.3.4)%*%X,t(X.4.4)%*%X,t(X.4.5)%*%X,t(X.4.6)%*%X),
                        cbind(t(X.1.5)%*%X,t(X.2.5)%*%X,t(X.3.5)%*%X,t(X.4.5)%*%X,t(X.5.5)%*%X,t(X.5.6)%*%X),
                        cbind(t(X.1.6)%*%X,t(X.2.6)%*%X,t(X.3.6)%*%X,t(X.4.6)%*%X,t(X.5.6)%*%X,t(X.6.6)%*%X))
      invV<-ginv(V)
    }
    ##Calculate the GMM loss.
    loss1<-as.vector(t(gbar)%*%invV%*%(gbar))      
    out1<-list("loss"=max(loss1*sum(n),loss1*sum(n)), "invV"=invV)
    out1
  }
  gmm.loss<-function(x,...) gmm.func(x,...)$loss
  
  ##Loss function for balance constraints, returns the squared imbalance along each dimension.
  bal.loss<-function(beta.curr){
    beta.curr<-matrix(beta.curr,k,no.treats-1)
    theta.curr<-X%*%beta.curr
    baseline.prob<-apply(theta.curr,1,function(x) (1+sum(exp(x)))^-1)
    probs.curr<-cbind(baseline.prob, exp(theta.curr[,1])*baseline.prob, exp(theta.curr[,2])*baseline.prob, exp(theta.curr[,3])*baseline.prob)
    probs.curr[,1]<-pmax(probs.min,probs.curr[,1])
    probs.curr[,2]<-pmax(probs.min,probs.curr[,2])
    probs.curr[,3]<-pmax(probs.min,probs.curr[,3])
    probs.curr[,4]<-pmax(probs.min,probs.curr[,4])
    norms<-apply(probs.curr,1,sum)
    probs.curr<-probs.curr/norms
    
    w.curr<-cbind(T1/probs.curr[,1] + T2/probs.curr[,2] - T3/probs.curr[,3] - T4/probs.curr[,4],
                  T1/probs.curr[,1] - T2/probs.curr[,2] - T3/probs.curr[,3] + T4/probs.curr[,4],
                  -T1/probs.curr[,1] + T2/probs.curr[,2] - T3/probs.curr[,3] + T4/probs.curr[,4])
    
    ##Generate mean imbalance.
    loss1<-sum(diag(abs(t(w.curr)%*%X%*%XprimeX.inv%*%t(X)%*%(w.curr))))
    loss1
  }
  
  gmm.gradient<-function(beta.curr, X.gmm=X, invV)
  {
    ##Designate a few objects in the function.
    beta.curr<-matrix(beta.curr,k,no.treats-1)
    X<-as.matrix(X.gmm)
    
    ##Designate sample size, number of treated and control observations,
    ##theta.curr, which are used to generate probabilities.
    ##Trim probabilities, and generate weights.
    n<-c(sum(T1), sum(T2), sum(T3), sum(T4))
    theta.curr<-X%*%beta.curr
    baseline.prob<-apply(theta.curr,1,function(x) (1+sum(exp(x)))^-1)
    probs.curr<-cbind(baseline.prob, exp(theta.curr[,1])*baseline.prob, exp(theta.curr[,2])*baseline.prob, exp(theta.curr[,3])*baseline.prob)
    probs.curr[,1]<-pmax(probs.min,probs.curr[,1])
    probs.curr[,2]<-pmax(probs.min,probs.curr[,2])
    probs.curr[,3]<-pmax(probs.min,probs.curr[,3])
    probs.curr[,4]<-pmax(probs.min,probs.curr[,4])
    norms<-apply(probs.curr,1,sum)
    probs.curr<-probs.curr/norms
    
    w.curr<-cbind(T1/probs.curr[,1] + T2/probs.curr[,2] - T3/probs.curr[,3] - T4/probs.curr[,4],
                  T1/probs.curr[,1] - T2/probs.curr[,2] - T3/probs.curr[,3] + T4/probs.curr[,4],
                  -T1/probs.curr[,1] + T2/probs.curr[,2] - T3/probs.curr[,3] + T4/probs.curr[,4])
    
    
    ##Generate the vector of mean imbalance by weights.
    w.curr.del<-1/sum(n)*t(X)%*%w.curr
    w.curr.del<-as.matrix(w.curr.del)
    w.curr<-as.matrix(w.curr)
    
    ##Generate g-bar, as in the paper.
    gbar<-c(1/sum(n)*t(X)%*%(T2-probs.curr[,2]), 
            1/sum(n)*t(X)%*%(T3-probs.curr[,3]),
            1/sum(n)*t(X)%*%(T4-probs.curr[,4]),
            w.curr.del)
    
    # TYPO?  I have a double negative in there.
    dgbar<-rbind(cbind(1/sum(n)*(t(apply(t(X),1,function(x) -x*probs.curr[,2]*(1-probs.curr[,2]))))%*%X,
                       1/sum(n)*(t(apply(t(X),1,function(x) x*probs.curr[,2]*probs.curr[,3])))%*%X,
                       1/sum(n)*(t(apply(t(X),1,function(x) x*probs.curr[,2]*probs.curr[,4])))%*%X,
                       1/sum(n)*(t(apply(t(X),1,function(x) x*(T1*probs.curr[,2]/probs.curr[,1] - T2*(1 - probs.curr[,2])/probs.curr[,2] - T3*probs.curr[,2]/probs.curr[,3] - T4*probs.curr[,2]/probs.curr[,4]))))%*%X,
                       1/sum(n)*(t(apply(t(X),1,function(x) x*(T1*probs.curr[,2]/probs.curr[,1] + T2*(1 - probs.curr[,2])/probs.curr[,2] - T3*probs.curr[,2]/probs.curr[,3] + T4*probs.curr[,2]/probs.curr[,4]))))%*%X,
                       1/sum(n)*(t(apply(t(X),1,function(x) x*(-T1*probs.curr[,2]/probs.curr[,1] - T2*(1 - probs.curr[,2])/probs.curr[,2] - T3*probs.curr[,2]/probs.curr[,3] + T4*probs.curr[,2]/probs.curr[,4]))))%*%X
    ),
    cbind(1/sum(n)*(t(apply(t(X),1,function(x) x*probs.curr[,2]*probs.curr[,3])))%*%X,
          1/sum(n)*(t(apply(t(X),1,function(x) -x*probs.curr[,3]*(1-probs.curr[,3]))))%*%X,
          1/sum(n)*(t(apply(t(X),1,function(x) x*probs.curr[,3]*probs.curr[,4])))%*%X,
          1/sum(n)*(t(apply(t(X),1,function(x) x*(T1*probs.curr[,3]/probs.curr[,1] + T2*probs.curr[,3]/probs.curr[,2] + T3*(1 - probs.curr[,3])/probs.curr[,3] - T4*probs.curr[,3]/probs.curr[,4]))))%*%X,
          1/sum(n)*(t(apply(t(X),1,function(x) x*(T1*probs.curr[,3]/probs.curr[,1] - T2*probs.curr[,3]/probs.curr[,2] + T3*(1 - probs.curr[,3])/probs.curr[,3] + T4*probs.curr[,3]/probs.curr[,4]))))%*%X,
          1/sum(n)*(t(apply(t(X),1,function(x) x*(-T1*probs.curr[,3]/probs.curr[,1] + T2*probs.curr[,3]/probs.curr[,2] + T3*(1 - probs.curr[,3])/probs.curr[,3] + T4*probs.curr[,3]/probs.curr[,4]))))%*%X
    ),
    cbind(1/sum(n)*(t(apply(t(X),1,function(x) x*probs.curr[,2]*probs.curr[,4])))%*%X,
          1/sum(n)*(t(apply(t(X),1,function(x) x*probs.curr[,3]*probs.curr[,4])))%*%X,
          1/sum(n)*(t(apply(t(X),1,function(x) -x*probs.curr[,4]*(1-probs.curr[,4]))))%*%X,
          1/sum(n)*(t(apply(t(X),1,function(x) x*(T1*probs.curr[,4]/probs.curr[,1] + T2*probs.curr[,4]/probs.curr[,2] - T3*probs.curr[,4]/probs.curr[,3] + T4*(1 - probs.curr[,4])/probs.curr[,4]))))%*%X,
          1/sum(n)*(t(apply(t(X),1,function(x) x*(T1*probs.curr[,4]/probs.curr[,1] - T2*probs.curr[,4]/probs.curr[,2] - T3*probs.curr[,4]/probs.curr[,3] - T4*(1 - probs.curr[,4])/probs.curr[,4]))))%*%X,
          1/sum(n)*(t(apply(t(X),1,function(x) x*(-T1*probs.curr[,4]/probs.curr[,1] + T2*probs.curr[,4]/probs.curr[,2] - T3*probs.curr[,4]/probs.curr[,3] - T4*(1 - probs.curr[,4])/probs.curr[,4]))))%*%X
    ))
    
    out<-sum(n)*2*dgbar%*%invV%*%gbar
    out
  }
  
  bal.gradient<-function(beta.curr)
  {
    beta.curr<-matrix(beta.curr,k,no.treats-1)
    theta.curr<-X%*%beta.curr
    baseline.prob<-apply(theta.curr,1,function(x) (1+sum(exp(x)))^-1)
    probs.curr<-cbind(baseline.prob, exp(theta.curr[,1])*baseline.prob, exp(theta.curr[,2])*baseline.prob, exp(theta.curr[,3])*baseline.prob)
    probs.curr[,1]<-pmax(probs.min,probs.curr[,1])
    probs.curr[,2]<-pmax(probs.min,probs.curr[,2])
    probs.curr[,3]<-pmax(probs.min,probs.curr[,3])
    probs.curr[,4]<-pmax(probs.min,probs.curr[,4])
    norms<-apply(probs.curr,1,sum)
    probs.curr<-probs.curr/norms
    
    w.curr<-cbind(T1/probs.curr[,1] + T2/probs.curr[,2] - T3/probs.curr[,3] - T4/probs.curr[,4],
                  T1/probs.curr[,1] - T2/probs.curr[,2] - T3/probs.curr[,3] + T4/probs.curr[,4],
                  -T1/probs.curr[,1] + T2/probs.curr[,2] - T3/probs.curr[,3] + T4/probs.curr[,4])
    
    dw.beta1<-cbind(t(apply(t(X),1,function(x) x*(T1*probs.curr[,2]/probs.curr[,1] - T2*(1 - probs.curr[,2])/probs.curr[,2] - T3*probs.curr[,2]/probs.curr[,3] - T4*probs.curr[,2]/probs.curr[,4]))),
                    t(apply(t(X),1,function(x) x*(T1*probs.curr[,2]/probs.curr[,1] + T2*(1 - probs.curr[,2])/probs.curr[,2] - T3*probs.curr[,2]/probs.curr[,3] + T4*probs.curr[,2]/probs.curr[,4]))),
                    t(apply(t(X),1,function(x) x*(-T1*probs.curr[,2]/probs.curr[,1] - T2*(1 - probs.curr[,2])/probs.curr[,2] - T3*probs.curr[,2]/probs.curr[,3] + T4*probs.curr[,2]/probs.curr[,4]))))
    
    dw.beta2<-cbind(t(apply(t(X),1,function(x) x*(T1*probs.curr[,3]/probs.curr[,1] + T2*probs.curr[,3]/probs.curr[,2] + T3*(1 - probs.curr[,3])/probs.curr[,3] - T4*probs.curr[,3]/probs.curr[,4]))),
                    t(apply(t(X),1,function(x) x*(T1*probs.curr[,3]/probs.curr[,1] - T2*probs.curr[,3]/probs.curr[,2] + T3*(1 - probs.curr[,3])/probs.curr[,3] + T4*probs.curr[,3]/probs.curr[,4]))),
                    t(apply(t(X),1,function(x) x*(-T1*probs.curr[,3]/probs.curr[,1] + T2*probs.curr[,3]/probs.curr[,2] + T3*(1 - probs.curr[,3])/probs.curr[,3] + T4*probs.curr[,3]/probs.curr[,4]))))
    
    dw.beta3<-cbind(t(apply(t(X),1,function(x) x*(T1*probs.curr[,4]/probs.curr[,1] + T2*probs.curr[,4]/probs.curr[,2] - T3*probs.curr[,4]/probs.curr[,3] + T4*(1 - probs.curr[,4])/probs.curr[,4]))),
                    t(apply(t(X),1,function(x) x*(T1*probs.curr[,4]/probs.curr[,1] - T2*probs.curr[,4]/probs.curr[,2] - T3*probs.curr[,4]/probs.curr[,3] - T4*(1 - probs.curr[,4])/probs.curr[,4]))),
                    t(apply(t(X),1,function(x) x*(-T1*probs.curr[,4]/probs.curr[,1] + T2*probs.curr[,4]/probs.curr[,2] - T3*probs.curr[,4]/probs.curr[,3] - T4*(1 - probs.curr[,4])/probs.curr[,4]))))
    
    loss1<-diag(t(w.curr)%*%X%*%XprimeX.inv%*%t(X)%*%(w.curr))
    
    out.1<-sapply(2*dw.beta1[,1:n]%*%X%*%XprimeX.inv%*%t(X)%*%(w.curr[,1]), 
                  function(x) ifelse((x > 0 & loss1[1] > 0) | (x < 0 & loss1[1] < 0),abs(x),-abs(x))) +
      sapply(2*dw.beta1[,(n+1):(2*n)]%*%X%*%XprimeX.inv%*%t(X)%*%(w.curr[,2]),
             function(x) ifelse((x > 0 & loss1[2] > 0) | (x < 0 & loss1[2] < 0),abs(x),-abs(x))) +
      sapply(2*dw.beta1[,(2*n+1):(3*n)]%*%X%*%XprimeX.inv%*%t(X)%*%(w.curr[,3]),
             function(x) ifelse((x > 0 & loss1[3] > 0) | (x < 0 & loss1[3] < 0),abs(x),-abs(x)))
    out.2<-sapply(2*dw.beta2[,1:n]%*%X%*%XprimeX.inv%*%t(X)%*%(w.curr[,1]), 
                  function(x) ifelse((x > 0 & loss1[1] > 0) | (x < 0 & loss1[1] < 0),abs(x),-abs(x))) +
      sapply(2*dw.beta2[,(n+1):(2*n)]%*%X%*%XprimeX.inv%*%t(X)%*%(w.curr[,2]),
             function(x) ifelse((x > 0 & loss1[2] > 0) | (x < 0 & loss1[2] < 0),abs(x),-abs(x))) +
      sapply(2*dw.beta2[,(2*n+1):(3*n)]%*%X%*%XprimeX.inv%*%t(X)%*%(w.curr[,3]),
             function(x) ifelse((x > 0 & loss1[3] > 0) | (x < 0 & loss1[3] < 0),abs(x),-abs(x)))
    out.3<-sapply(2*dw.beta3[,1:n]%*%X%*%XprimeX.inv%*%t(X)%*%(w.curr[,1]), 
                  function(x) ifelse((x > 0 & loss1[1] > 0) | (x < 0 & loss1[1] < 0),abs(x),-abs(x))) +
      sapply(2*dw.beta3[,(n+1):(2*n)]%*%X%*%XprimeX.inv%*%t(X)%*%(w.curr[,2]),
             function(x) ifelse((x > 0 & loss1[2] > 0) | (x < 0 & loss1[2] < 0),abs(x),-abs(x))) +
      sapply(2*dw.beta3[,(2*n+1):(3*n)]%*%X%*%XprimeX.inv%*%t(X)%*%(w.curr[,3]),
             function(x) ifelse((x > 0 & loss1[3] > 0) | (x < 0 & loss1[3] < 0),abs(x),-abs(x)))
    out<-c(out.1, out.2, out.3)
    out
  }
  
  n<-length(treat)
  n.t<-c(sum(T1),sum(T2),sum(T3),sum(T4))
  x.orig<-x<-cbind(as.matrix(X))
  
  ##Run multionmial logit
  dat.dummy<-data.frame(treat=treat,X)
  #Need to generalize for different dimensioned X's
  xnam<- paste0("X", 1:dim(X)[2])
  fmla <- as.formula(paste("as.factor(treat) ~ -1 + ", paste(xnam, collapse= "+")))
  mnl1<-multinom(fmla, data=dat.dummy,trace=FALSE)
  mcoef<-matrix(coef(mnl1),k,no.treats-1,byrow=TRUE)
  mcoef[is.na(mcoef[,1]),1]<-0
  mcoef[is.na(mcoef[,2]),2]<-0
  mcoef[is.na(mcoef[,3]),3]<-0
  probs.mnl<-cbind(1/(1+exp(X%*%mcoef[,1])+exp(X%*%mcoef[,2])+exp(X%*%mcoef[,3])),
                   exp(X%*%mcoef[,1])/(1+exp(X%*%mcoef[,1])+exp(X%*%mcoef[,2])+exp(X%*%mcoef[,3])),
                   exp(X%*%mcoef[,2])/(1+exp(X%*%mcoef[,1])+exp(X%*%mcoef[,2])+exp(X%*%mcoef[,3])),
                   exp(X%*%mcoef[,3])/(1+exp(X%*%mcoef[,1])+exp(X%*%mcoef[,2])+exp(X%*%mcoef[,3])))
  colnames(probs.mnl)<-c("p1","p2","p3","p4")
  probs.mnl[,1]<-pmax(probs.min,probs.mnl[,1])
  probs.mnl[,2]<-pmax(probs.min,probs.mnl[,2])
  probs.mnl[,3]<-pmax(probs.min,probs.mnl[,3])
  probs.mnl[,4]<-pmax(probs.min,probs.mnl[,4])
  norms<-apply(probs.mnl,1,sum)
  probs.mnl<-probs.mnl/norms
  mnl1$fit<-matrix(probs.mnl,nrow=n,ncol=no.treats)
  beta.curr<-mcoef
  beta.curr[is.na(beta.curr)]<-0
  
  
  alpha.func<-function(alpha) gmm.loss(beta.curr*alpha)
  beta.curr<-beta.curr*optimize(alpha.func,interval=c(.8,1.1))$min
  
  ##Generate estimates for balance and CBPSE
  gmm.init<-beta.curr
  
  if(twostep)
  {
    this.invV<-gmm.func(gmm.init)$invV
  }
  
  if (twostep)
  {
    opt.bal<-optim(gmm.init, bal.loss, control=list("maxit"=iterations), method="BFGS", gr = bal.gradient, hessian=TRUE)
  }
  else
  {
    opt.bal<-tryCatch({
      optim(gmm.init, bal.loss, control=list("maxit"=iterations), method="BFGS", hessian=TRUE)
    },
    error = function(err)
    {
      return(optim(gmm.init, bal.loss, control=list("maxit"=iterations), method="Nelder-Mead", hessian=TRUE))
    })	  
  }
  beta.bal<-opt.bal$par	  
  if(bal.only) opt1<-opt.bal 
  
  if(twostep)
  {
    if (gmm.loss(gmm.init) < gmm.loss(beta.bal))
    {
      this.invV<-gmm.func(gmm.init)$invV
    }
    else
    {
      this.invV<-gmm.func(beta.bal)$invV
    }
    if(bal.only)
    {
      this.invV<-gmm.func(beta.bal)$invV
    }
  }
  
  if(!bal.only)
  {
    if (twostep)
    {			        
      gmm.glm.init<-optim(gmm.init, gmm.loss, control=list("maxit"=iterations), method="BFGS", gr = gmm.gradient, hessian=TRUE, invV = this.invV)
      gmm.bal.init<-optim(beta.bal, gmm.loss, control=list("maxit"=iterations), method="BFGS", gr = gmm.gradient, hessian=TRUE, invV = this.invV)
    }
    else
    {
      gmm.glm.init<-tryCatch({
        optim(gmm.init,gmm.loss, control=list("maxit"=iterations), method="BFGS", hessian=TRUE)
      },
      error = function(err)
      {
        return(optim(gmm.init,gmm.loss, control=list("maxit"=iterations), method="Nelder-Mead", hessian=TRUE))
      })
      gmm.bal.init<-tryCatch({
        optim(beta.bal, gmm.loss, control=list("maxit"=iterations), method="BFGS", hessian=TRUE)
      },
      error = function(err)
      {
        return(optim(beta.bal, gmm.loss, control=list("maxit"=iterations), method="Nelder-Mead", hessian=TRUE))
      })
    }
    
    if(gmm.glm.init$val<gmm.bal.init$val) opt1<-gmm.glm.init else opt1<-gmm.bal.init
  }
  
  ##Generate probabilities
  beta.opt<-matrix(opt1$par,nrow=k,ncol=no.treats-1)
  theta.opt<-X%*%beta.opt
  baseline.prob<-apply(theta.opt,1,function(x) (1+sum(exp(x)))^-1)
  probs.opt<-cbind(baseline.prob, exp(theta.opt[,1])*baseline.prob, exp(theta.opt[,2])*baseline.prob, exp(theta.opt[,3])*baseline.prob)
  probs.opt[,1]<-pmax(probs.min,probs.opt[,1])
  probs.opt[,2]<-pmax(probs.min,probs.opt[,2])
  probs.opt[,3]<-pmax(probs.min,probs.opt[,3])
  probs.opt[,4]<-pmax(probs.min,probs.opt[,4])
  norms<-apply(probs.opt,1,sum)
  probs.opt<-probs.opt/norms
  
  J.opt<-ifelse(twostep, gmm.func(beta.opt, invV = this.invV)$loss, gmm.func(beta.opt)$loss)
  
  if ((J.opt > gmm.loss(mcoef)) & (bal.loss(beta.opt) > bal.loss(mcoef)))
  {	  
    beta.opt<-mcoef
    probs.opt<-probs.mnl
    J.opt <- gmm.loss(mcoef)
    warning("Optimization failed.  Results returned are for MLE.")
  }
  
  ##Generate weights
  w.opt<-cbind(T1/probs.opt[,1] + T2/probs.opt[,2] - T3/probs.opt[,3] - T4/probs.opt[,4],
               T1/probs.opt[,1] - T2/probs.opt[,2] - T3/probs.opt[,3] + T4/probs.opt[,4],
               -T1/probs.opt[,1] + T2/probs.opt[,2] - T3/probs.opt[,3] + T4/probs.opt[,4])
  
  #How are residuals now defined?
  residuals<-cbind(T1-probs.opt[,1],T2-probs.opt[,2],T3-probs.opt[,3],T4-probs.opt[,4])
  deviance <- -2*c(sum(T1*log(probs.opt[,1])+T2*log(probs.opt[,2])+T3*log(probs.opt[,3])+T4*log(probs.opt[,4])))
  nulldeviance <- -2*c(sum(T1*log(mean(T1))+T2*log(mean(T2))+T3*log(mean(T3))+T4*log(mean(T4))))
  
  norm1<-norm2<-norm3<-norm4<-1
  if (standardize)
  {
    norm1<-sum(1/probs.opt[which(T1==1),1])
    norm2<-sum(1/probs.opt[which(T2==1),2])
    norm3<-sum(1/probs.opt[which(T3==1),3])
    norm4<-sum(1/probs.opt[which(T4==1),4])
  }
  
  w<-(T1 == 1)/probs.opt[,1]/norm1 + (T2 == 1)/probs.opt[,2]/norm2 + (T3 == 1)/probs.opt[,3]/norm3 + (T4 == 1)/probs.opt[,4]/norm4
  
  if (twostep)
  {
    W<-this.invV
  }
  else
  {
    W<-gmm.func(beta.opt)$invV
  }
  X.G.1.1<-t(-X*probs.opt[,2]*(1-probs.opt[,2]))%*%X
  X.G.1.2<-t(X*probs.opt[,2]*probs.opt[,3])%*%X
  X.G.1.3<-t(X*probs.opt[,2]*probs.opt[,4])%*%X
  X.G.1.4<-t(X*probs.opt[,2]*(T1/probs.opt[,1] - T2*(1-probs.opt[,2])/probs.opt[,2]^2 - T3/probs.opt[,3] - T4/probs.opt[,4]))%*%X
  X.G.1.5<-t(X*probs.opt[,2]*(T1/probs.opt[,1] + T2*(1-probs.opt[,2])/probs.opt[,2]^2 - T3/probs.opt[,3] + T4/probs.opt[,4]))%*%X
  X.G.1.6<-t(X*probs.opt[,2]*(-T1/probs.opt[,1] - T2*(1-probs.opt[,2])/probs.opt[,2]^2 - T3/probs.opt[,3] + T4/probs.opt[,4]))%*%X
  X.G.2.1<-t(X*probs.opt[,2]*probs.opt[,3])%*%X
  X.G.2.2<-t(-X*probs.opt[,3]*(1-probs.opt[,3]))%*%X
  X.G.2.3<-t(X*probs.opt[,3]*probs.opt[,4])%*%X
  X.G.2.4<-t(X*probs.opt[,3]*(T1/probs.opt[,1] + T2/probs.opt[,2] + T3*(1-probs.opt[,3])/probs.opt[,3]^2 - T4/probs.opt[,4]))%*%X
  X.G.2.5<-t(X*probs.opt[,3]*(T1/probs.opt[,1] - T2/probs.opt[,2] + T3*(1-probs.opt[,3])/probs.opt[,3]^2 + T4/probs.opt[,4]))%*%X
  X.G.2.6<-t(X*probs.opt[,3]*(-T1/probs.opt[,1] + T2/probs.opt[,2] + T3*(1-probs.opt[,3])/probs.opt[,3]^2 + T4/probs.opt[,4]))%*%X
  X.G.3.1<-t(X*probs.opt[,2]*probs.opt[,4])%*%X
  X.G.3.2<-t(X*probs.opt[,3]*probs.opt[,4])%*%X
  X.G.3.3<-t(-X*probs.opt[,4]*(1-probs.opt[,4]))%*%X
  X.G.3.4<-t(X*probs.opt[,4]*(T1/probs.opt[,1] + T2/probs.opt[,2] - T3/probs.opt[,3] + T4*(1-probs.opt[,4])/probs.opt[,4]^2))%*%X
  X.G.3.5<-t(X*probs.opt[,4]*(T1/probs.opt[,1] - T2/probs.opt[,2] - T3/probs.opt[,3] - T4*(1-probs.opt[,4])/probs.opt[,4]^2))%*%X
  X.G.3.6<-t(X*probs.opt[,4]*(-T1/probs.opt[,1] + T2/probs.opt[,2] - T3/probs.opt[,3] - T4*(1-probs.opt[,4])/probs.opt[,4]^2))%*%X
  XW.1<-X*(T2-probs.opt[,2])
  XW.2<-X*(T3-probs.opt[,3])
  XW.3<-X*(T4-probs.opt[,4])
  XW.4<-X*( T1/probs.opt[,1] + T2/probs.opt[,2] - T3/probs.opt[,3] - T4/probs.opt[,4])
  XW.5<-X*( T1/probs.opt[,1] - T2/probs.opt[,2] - T3/probs.opt[,3] + T4/probs.opt[,4])
  XW.6<-X*(-T1/probs.opt[,1] + T2/probs.opt[,2] - T3/probs.opt[,3] + T4/probs.opt[,4])
  G<-1/n*rbind(cbind(X.G.1.1,X.G.1.2,X.G.1.3,X.G.1.4,X.G.1.5,X.G.1.6),
               cbind(X.G.2.1,X.G.2.2,X.G.2.3,X.G.2.4,X.G.2.5,X.G.2.6),
               cbind(X.G.3.1,X.G.3.2,X.G.3.3,X.G.3.4,X.G.3.5,X.G.3.6))
  W1<-rbind(t(XW.1),t(XW.2),t(XW.3),t(XW.4),t(XW.5),t(XW.6))
  Omega<-1/n*(W1%*%t(W1))
  vcov<-ginv(G%*%W%*%t(G))%*%G%*%W%*%Omega%*%W%*%t(G)%*%ginv(G%*%W%*%t(G))
  
  colnames(probs.opt)<-treat.names
  
  class(beta.opt) <- "coef"
  
  output<-list("coefficients"=beta.opt,"fitted.values"=probs.opt,"deviance"=deviance,"weights"=w, 
               "y"=treat,"x"=X,"converged"=opt1$conv,"J"=J.opt,"var"=vcov,   
               "mle.J"=ifelse(twostep, gmm.func(mcoef, invV = this.invV)$loss, gmm.loss(mcoef)))
  
  class(output)<- c("CBPS")
  output
}

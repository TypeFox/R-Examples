CBPS.Continuous<-function(treat, X, X.bal, method, k, XprimeX.inv, bal.only, iterations, standardize, twostep, ...)
{
	probs.min<-1e-6
	
	##The gmm objective function--given a guess of beta, constructs the GMM J statistic.
	  gmm.func<-function(params.curr,X.gmm=Xtilde,invV=NULL){
		##Generate probabilities.
		##Trim probabilities, and generate weights.
		beta.curr<-params.curr[-length(params.curr)]
		sigmasq<-exp(params.curr[length(params.curr)])
		
		probs.curr<-dnorm(Ttilde, mean = Xtilde%*%beta.curr, sd = sqrt(sigmasq), log = TRUE)
		probs.curr<-pmin(log(1-probs.min),probs.curr)
		probs.curr<-pmax(log(probs.min),probs.curr)
		
		w.curr<-Ttilde*exp(stabilizers - probs.curr)
		
		##Generate the vector of mean imbalance by weights.
		w.curr.del<-1/n*t(Xtilde)%*%w.curr
		w.curr.del<-as.matrix(w.curr.del)
		w.curr<-as.matrix(w.curr)

		##Generate g-bar, as in the paper.
		gbar<-c(1/n*t(Xtilde)%*%(Ttilde-Xtilde%*%beta.curr)/sigmasq,
				1/n*t(n.identity.vec)%*%((Ttilde - Xtilde%*%beta.curr)^2/(2*sigmasq^2) - 1/(2*sigmasq)),
				w.curr.del)
				
		##Generate the covariance matrix used in the GMM estimate.
		if (is.null(invV))
		{
			Xtilde.1.1<-1/sigmasq*t(Xtilde)%*%(Xtilde)
			Xtilde.1.2<-0*t(Xtilde)%*%n.identity.vec
			Xtilde.1.3<-t(Xtilde)%*%(Xtilde)
			Xtilde.2.2<-t(n.identity.vec)%*%n.identity.vec*(2*sigmasq^2)^(-1)
			Xtilde.2.3<-t(Xtilde)%*%(-Xtilde%*%beta.curr/sigmasq)
			Xtilde.3.3<-t(Xtilde)%*%sweep(Xtilde,MARGIN=1,pmin(as.vector(exp((Xtilde%*%beta.curr)^2/sigmasq + log(sigmasq + (Xtilde%*%beta.curr)^2))), 10^250),'*')
						
			V<-rbind(1/n*cbind(Xtilde.1.1,Xtilde.1.2,Xtilde.1.3),
					 1/n*cbind(t(Xtilde.1.2),Xtilde.2.2,t(Xtilde.2.3)),
					 1/n*cbind(Xtilde.1.3,Xtilde.2.3,Xtilde.3.3))
			invV<-ginv(V)
		}
		
		##Calculate the GMM loss.
		loss1<-t(gbar)%*%invV%*%(gbar)
		out1<-list("loss"=loss1*n, "invV"=invV)
		out1
	  }
	  gmm.loss<-function(x,...) gmm.func(x,...)$loss


	  
	  ##Loss function for balance constraints, returns the squared imbalance along each dimension.
	  bal.loss<-function(params.curr){
	  	beta.curr<-params.curr[-length(params.curr)]
		  sigmasq<-exp(params.curr[length(params.curr)])
		
		  probs.curr<-dnorm(Ttilde, mean = Xtilde%*%beta.curr, sd = sqrt(sigmasq), log = TRUE)
		  probs.curr<-pmin(log(1-probs.min),probs.curr)
		  probs.curr<-pmax(log(probs.min),probs.curr)
		
		  w.curr<-Ttilde*exp(stabilizers - probs.curr)
		
		  ##Generate mean imbalance.
		  loss1<-t(w.curr)%*%Xtilde%*%XtildeprimeXtilde.inv%*%t(Xtilde)%*%(w.curr)
		  loss1
	  }
	  
	  gmm.gradient<-function(params.curr, X.gmm=X, invV)
	  {
	  	##Generate probabilities.
		##Trim probabilities, and generate weights.
		beta.curr<-params.curr[-length(params.curr)]
		sigmasq<-exp(params.curr[length(params.curr)])
			
		probs.curr<-dnorm(Ttilde, mean = Xtilde%*%beta.curr, sd = sqrt(sigmasq), log = TRUE)
		probs.curr<-pmin(log(1-probs.min),probs.curr)
		probs.curr<-pmax(log(probs.min),probs.curr)
		
		w.curr<-Ttilde*exp(stabilizers - probs.curr)
		

		##Generate the vector of mean imbalance by weights.
		w.curr.del<-1/n*t(Xtilde)%*%w.curr
		w.curr.del<-as.matrix(w.curr.del)
		w.curr<-as.matrix(w.curr)

		##Generate g-bar, as in the paper.
		gbar<-c(1/n*t(Xtilde)%*%(Ttilde-Xtilde%*%beta.curr)/sigmasq,
				1/n*t(n.identity.vec)%*%((Ttilde - Xtilde%*%beta.curr)^2/(2*sigmasq^2) - 1/(2*sigmasq)),
				w.curr.del)
		
		dgbar.1.1<-t(-Xtilde)%*%Xtilde/sigmasq
		dgbar.1.2<-matrix(-(Ttilde - Xtilde%*%beta.curr)/(sigmasq^2), nrow = 1)%*%Xtilde
		dgbar.2.2<-t(n.identity.vec)%*%(-(Ttilde - Xtilde%*%beta.curr)^2/(sigmasq^3) + (2*sigmasq^2)^(-1))
		dgbar.3.1<-sweep(t(Xtilde), MARGIN=2, -(Ttilde-Xtilde%*%beta.curr)/sigmasq*w.curr,'*')%*%Xtilde
		dgbar.3.2<-matrix(-w.curr*((Ttilde - Xtilde%*%beta.curr)^2/(2*sigmasq^2) - 1/(2*sigmasq)), nrow = 1)%*%Xtilde
		dgbar<-1/n*cbind(rbind(dgbar.1.1, dgbar.1.2*sigmasq),
						 rbind(t(dgbar.1.2), dgbar.2.2*sigmasq),
						 rbind(dgbar.3.1, dgbar.3.2*sigmasq))
		
		out<-2*n*dgbar%*%invV%*%gbar
		out
	  }
	  
	  bal.gradient<-function(params.curr, invV=NULL)
	  {
	  	beta.curr<-params.curr[-length(params.curr)]
		  sigmasq<-exp(params.curr[length(params.curr)])
		
	  	probs.curr<-dnorm(Ttilde, mean = Xtilde%*%beta.curr, sd = sqrt(sigmasq), log = TRUE)
		  probs.curr<-pmin(log(1-probs.min),probs.curr)
		  probs.curr<-pmax(log(probs.min),probs.curr)
		
		  w.curr<-Ttilde*exp(stabilizers - probs.curr)
		
		  dw<-1/n*rbind(sweep(t(Xtilde), MARGIN=2, -(Ttilde-Xtilde%*%beta.curr)/sigmasq*w.curr,'*'),
			         		  sigmasq*as.vector(-w.curr*((Ttilde - Xtilde%*%beta.curr)^2/(2*sigmasq^2) - 1/(2*sigmasq))))
		  ##Generate mean imbalance.
      out<-2*n*dw%*%Xtilde%*%XtildeprimeXtilde.inv%*%t(Xtilde)%*%(w.curr)
      out
	  }
	  
	  n<-length(treat)
	  x.orig<-x<-cbind(as.matrix(X))
	  Xtilde<-apply(X,2,function(x) (x - mean(x)))
	  XtildeprimeXtilde.inv<-ginv(t(Xtilde)%*%Xtilde)
	  Ttilde<-(treat-mean(treat))
	  n.identity.vec<-matrix(1,nrow=n,ncol=1)
  
	##Run linear regression
	  lm1<-lm(Ttilde ~ -1 + Xtilde)
	  mcoef<-coef(lm1)
	  mcoef[is.na(mcoef)]<-0
	  sigmasq<-mean((Ttilde - Xtilde%*%mcoef)^2)
	  Ttilde.hat<-apply(Xtilde, 1, function(x) x%*%mcoef)
	  stabilizers<-log(sapply(Ttilde, function(t) mean(pmin(pmax(dnorm(t, mean = Ttilde.hat, sd = sqrt(sigmasq)), probs.min), 1-probs.min))))
	  
	  probs.mle<-dnorm(Ttilde, mean = Xtilde%*%mcoef, sd = sqrt(sigmasq), log = TRUE)
	  probs.mle<-pmin(log(1-probs.min),probs.mle)
	  probs.mle<-pmax(log(probs.min),probs.mle)
	  
	  params.curr<-c(mcoef, log(sigmasq))
	  mle.J<-gmm.loss(params.curr)
	  mle.bal<-bal.loss(params.curr)

    alpha.func<-function(alpha) gmm.loss(params.curr*alpha)      
	  params.curr<-params.curr*optimize(alpha.func,interval=c(.8,1.1))$min
  
	  if (twostep)
	  {
	    glm.invV<-gmm.func(params.curr)$invV
  	}
	  
	  ##Generate estimates for balance and CBPS
	  gmm.init<-params.curr
	
	  if (twostep)
	  {
		opt.bal<-optim(gmm.init, bal.loss, control=list("maxit"=iterations), method="BFGS", gr = bal.gradient, hessian = TRUE)
	  }
	  else
	  {
		opt.bal<-tryCatch({optim(gmm.init, bal.loss, control=list("maxit"=iterations), method="BFGS", hessian=TRUE)},
						   error = function(err) {return(optim(gmm.init, bal.loss, control=list("maxit"=iterations), method="Nelder-Mead", hessian=TRUE))})
	  }
	  params.bal<-opt.bal$par
		
	  if(bal.only) opt1<-opt.bal
	  
    pick.glm<-0
	  if(!bal.only)
	  {
		if (twostep)
		{
			gmm.glm.init<-optim(gmm.init, gmm.loss, control=list("maxit"=iterations), method="BFGS", hessian=TRUE, invV = glm.invV, gr = gmm.gradient)
			gmm.bal.init<-optim(params.bal, gmm.loss, control=list("maxit"=iterations), method="BFGS", hessian=TRUE, invV = glm.invV, gr = gmm.gradient)
		}
		else	
		{
			gmm.glm.init<-tryCatch({optim(params.curr, gmm.loss, control=list("maxit"=iterations), method="BFGS", hessian=TRUE)},
									error = function(err) {return(optim(params.curr, gmm.loss, control=list("maxit"=iterations), method="Nelder-Mead", hessian=TRUE))})
			gmm.bal.init<-tryCatch({optim(params.bal, gmm.loss, control=list("maxit"=iterations), method="BFGS", hessian=TRUE)},
								error = function(err) {return(optim(params.bal, gmm.loss, control=list("maxit"=iterations), method="Nelder-Mead", hessian=TRUE))})
		}	
		if(gmm.glm.init$val<gmm.bal.init$val) opt1<-gmm.glm.init else opt1<-gmm.bal.init
		pick.glm<-ifelse(gmm.glm.init$val<gmm.bal.init$val,1,0)
	  }
	  				
	  ##Generate probabilities
	  params.opt<-opt1$par
	  beta.opt<-params.opt[-length(params.opt)]
	  sigmasq<-exp(params.opt[length(params.opt)])
	  probs.opt<-dnorm(Ttilde, mean = Xtilde%*%beta.opt, sd = sqrt(sigmasq), log = TRUE)
	  probs.opt<-pmin(log(1-probs.min),probs.opt)
	  probs.opt<-pmax(log(probs.min),probs.opt)
  
	  J.opt<-ifelse(twostep, gmm.func(params.opt, invV = glm.invV)$loss, gmm.func(params.opt)$loss)

	  if ((J.opt > mle.J) & (bal.loss(params.opt) > mle.bal))
		{	  
			beta.opt<-mcoef
			probs.opt<-probs.mle
			J.opt <-mle.J
			warning("Optimization failed.  Results returned are for MLE.")
		}
		
	  ##Generate weights
	  w.opt<-exp(stabilizers - probs.opt)
    if(standardize) w.opt<-w.opt/sum(w.opt)
		
	  #How are residuals now defined?
	  residuals<- Ttilde - Xtilde%*%beta.opt
	  deviance <- -2*sum(probs.opt)
	  
	  if (twostep)
	  {
		W<-glm.invV
	  }
	  else
	  {
		W<-gmm.func(params.opt)$invV
	  }
	  
	  XG.1.1<-t(-Xtilde)%*%Xtilde/sigmasq
	  XG.1.2<-t(-Xtilde)%*%(Ttilde - Xtilde%*%beta.opt)/(sigmasq^2)
	  XG.1.3<-sweep(t(Xtilde), MARGIN=2, -(Ttilde-Xtilde%*%beta.opt)/sigmasq*Ttilde*w.opt,'*')%*%Xtilde
	  XG.2.2<-t(n.identity.vec)%*%(-(Ttilde - Xtilde%*%beta.opt)^2/(sigmasq^3) + (2*sigmasq^2)^(-1))
	  XG.2.3<-matrix(-Ttilde*w.opt*((Ttilde - Xtilde%*%beta.opt)^2/(2*sigmasq^2) - 1/(2*sigmasq)), nrow = 1)%*%Xtilde
	  
	  XW.1<-sweep(Xtilde,MARGIN=1,as.vector((Ttilde-Xtilde%*%beta.opt)/sigmasq),'*')
	  XW.2<-as.vector((Ttilde - Xtilde%*%beta.opt)^2/(2*sigmasq^2) - 1/(2*sigmasq))
	  XW.3<-sweep(Xtilde,MARGIN=1,as.vector(Ttilde*exp(stabilizers - probs.opt)),'*')
	  	  
	  G<-1/n*rbind(cbind(XG.1.1,XG.1.2,XG.1.3),
				   cbind(t(XG.1.2),XG.2.2,XG.2.3))
	  W1<-rbind(t(XW.1),t(XW.2),t(XW.3))
	  Omega<-1/n*(W1%*%t(W1))
	  
	  vcov<-(ginv(G%*%W%*%t(G))%*%G%*%W%*%Omega%*%W%*%t(G)%*%ginv(G%*%W%*%t(G)))[1:k,1:k] 
	  #Reverse the centering from using Ttilde and Xtilde
	  beta.tilde<-beta.opt
	  beta.opt<-ginv(t(X)%*%X)%*%t(X)%*%(Xtilde%*%beta.opt + mean(treat))
	  vcov<-vcov + ginv(t(Xtilde)%*%Xtilde)%*%(t(Xtilde)%*%Xtilde)*as.numeric(2*t(X%*%beta.opt)%*%n.identity.vec*mean(treat) - mean(treat)^2 + t(Xtilde%*%beta.tilde)%*%(Xtilde%*%beta.tilde) - t(X%*%beta.opt)%*%(X%*%beta.opt))
	  
	  class(beta.opt)<-"coef"
		
	  output<-list("coefficients"=beta.opt, "sigmasq"=sigmasq,"fitted.values"=dnorm(treat,X%*%beta.opt,sigmasq),"deviance"=deviance,
                 "weights"=w.opt,"y"=treat,"x"=X,"converged"=opt1$conv,"J"=J.opt,"var"=vcov, "mle.J"=mle.J)
  
	  class(output)<- c("CBPSContinuous","CBPS")
	  output
}

####


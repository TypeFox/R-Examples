mle_marginal <-
function(x,y,R,S,family,exposure=rep(1,length(y)),sd.error=FALSE,zt=TRUE){
    # fit gamma-model
    n<-length(x)
    my.gamma<-glm(x~-1 +R,family=Gamma(link="log"))
    sd.alpha=NULL
    if (sd.error==TRUE){
		sd.alpha<-sqrt(diag(vcov(my.gamma)))
	}
    alpha<-my.gamma$coefficients
    delta<-summary(my.gamma)$dispersion
    mu<-as.vector(exp(R%*%alpha))
	# fit ztp/poisson model
    if (zt==TRUE){
        ztp.model<-ztp.glm(y,S,exposure,sd.error)
	   beta<-ztp.model$coefficients
	   sd.beta<-ztp.model$sd
    }
    if (zt==FALSE){
        pois.model=glm(y~S-1,offset=log(exposure),family=poisson(link="log"))
        beta<-coef(pois.model)
        sd.beta=NULL
        if (sd.error==TRUE){
            sd.beta<-sqrt(diag(vcov(pois.model)))
        }
    }
    lambda<-as.vector(exp(S%*%beta))*exposure
    # inference by margins
	theta_initial<-BiCopEst(rank(x-mu)/(length(x)+1),rank(y-lambda)/(length(y)+1),family=family)$par
	tau_initial=BiCopPar2Tau(par=theta_initial,family=family)
	u<-pgam(x,mu,delta)
    if (zt==TRUE){
	   v<-pztp(y,lambda)
	   vv<-pztp(y-1,lambda)
    }
    if (zt==FALSE){
        v<-ppois(y,lambda)
	   vv<-ppois(y-1,lambda)
    }
	foo<-function(para){
		theta0<-z2theta(para,family)
		out<-(-sum(log(D_u(u,v,theta0,family)- D_u(u,vv,theta0,family))))
		return(out)
	}
	para_initial<-theta2z(theta_initial,family)
	para.ifm<-optim(para_initial,foo,method="BFGS")$par
	theta.ifm<-z2theta(para.ifm,family)
	tau.ifm<-BiCopPar2Tau(par=theta.ifm,family=family)
    	theta=theta_initial
    	tau=BiCopPar2Tau(par=theta,family=family)
    loglik<-loglik_joint(alpha,beta,0,delta,x,y,R,S,1,exposure,zt)
	loglik.ifm<-loglik_joint(alpha,beta,theta.ifm,delta,x,y,R,S,family,exposure,zt)
mu<-as.vector(exp(R%*%alpha))
    lambda<-as.vector(exp(S%*%beta)*exposure)
    ll<-log(density_joint(x,y,mu,delta,lambda,0,1,zt))
    ll.ifm<-log(density_joint(x,y,mu,delta,lambda,theta.ifm,family,zt))
    return(list(alpha=alpha,beta=beta,sd.alpha=sd.alpha,sd.beta=sd.beta,delta=delta,theta=0,family0=family,theta.ifm=theta.ifm,tau.ifm=tau.ifm,loglik=loglik,loglik.ifm=loglik.ifm,family=1,ll=ll,ll.ifm=ll.ifm))
    
}

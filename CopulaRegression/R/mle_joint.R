mle_joint <-
function(alpha0,beta0,theta0,delta0,x,y,R,S,family,exposure=rep(1,length(y)),sd.error=FALSE,zt=TRUE){
   p<-ncol(R)
   q<-ncol(S)
    foo<-function(para,x,y,R,S,family,exposure,zt){
	 p<-ncol(R)
    	q<-ncol(S)
        alpha<-para[1:p]
	beta<-para[(p+1):(p+q)]
	theta<-z2theta(para[p+q+1],family)
	#theta<-para[p+q+1]
    delta<-exp(para[p+q+2])
	return(loglik_joint(alpha,beta,theta,delta,x,y,R,S,family,exposure,zt=zt))
	}
	para0<-c(alpha0,beta0,theta2z(theta0,family),log(delta0))	
	dummy<-optim(para0,foo,x=x,y=y,R=R,S=S,family=family,exposure=exposure,zt=zt,method="BFGS",hessian=sd.error)
	
    # transform parameters back
    out<-dummy$par
    alpha<-out[1:p]
    beta<-out[(p+1):(p+q)]
    theta=z2theta(out[p+q+1],family)
	#theta<-out[p+q+1]
    delta<-exp(out[p+q+2])
	sd.alpha=NULL
	sd.beta=NULL
    sd.g.theta=NULL
    if (sd.error==TRUE){
    	hessian<-dummy$hessian
   	Hinv<-ginv(hessian)
	sd<-sqrt(diag(Hinv))
	sd.alpha<-sd[1:p]
	sd.beta<-sd[(p+1):(p+q)]
	sd.g.theta<-sd[(p+q+1)]
	}
    tau<-BiCopPar2Tau(par=theta,family=family)
	loglik<-loglik_joint(alpha,beta,theta,delta,x,y,R,S,family,zt)
    mu<-as.vector(exp(R%*%alpha))
    lambda<-as.vector(exp(S%*%beta)*exposure)
    ll<-log(density_joint(x,y,mu,delta,lambda,theta,family,zt))
return(list(alpha=alpha,beta=beta,sd.alpha=sd.alpha,sd.beta=sd.beta,sd.g.theta=sd.g.theta,delta=delta,tau=tau,theta=theta,family=family,loglik=loglik,ll=ll))
}

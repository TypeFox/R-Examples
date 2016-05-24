test.equality <- function(y, x = NULL, arbmean=TRUE, arbvar=FALSE, mu=NULL, sigma=NULL, beta=NULL, lambda=NULL,...){

if(arbmean==arbvar) stop("Change either 'arbmean' or 'arbvar'!")

if(arbmean==FALSE){
	w=1
	while(w==1){
	if(is.null(x)){
	H0=normalmixEM(x=y,arbmean=FALSE,arbvar=TRUE, mu=mu, sigma=sigma, lambda=lambda,...)
	k=length(H0$lambda)
#	H1=normalmixEM(x=y,arbmean=TRUE,arbvar=TRUE,lambda=H0$lambda,mu=rep(H0$mu,k)*(1:k),sigma=(H0$scale*H0$sigma),...)
	H1=normalmixEM(x=y,arbmean=TRUE,arbvar=TRUE,lambda=H0$lambda,mu=NULL,sigma=(H0$scale*H0$sigma),...)
	D=2*(H1$loglik-H0$loglik)
	df=k-1
	alpha=1-pchisq(D,df=df)
	} else{
	H0=regmixEM(y=y,x=x,arbmean=FALSE,arbvar=TRUE,beta=beta, sigma=sigma, lambda=lambda,...)
	k=length(H0$lambda)
#	H1=regmixEM(y=y,x=x,arbmean=TRUE,arbvar=TRUE,lambda=H0$lambda,beta=matrix(rep(H0$beta,k),k)*(1:k),sigma=(H0$scale*H0$sigma),...)
	H1=regmixEM(y=y,x=x,arbmean=TRUE,arbvar=TRUE,lambda=H0$lambda,beta=NULL,sigma=(H0$scale*H0$sigma),...)
	p=nrow(H1$beta)
	D=2*(H1$loglik-H0$loglik)
	df=p*(k-1)
	alpha=1-pchisq(D,df=df)
	}
	if(D<0){
	w=1 
	mu=NULL
	sigma=NULL
	lambda=NULL
	} else w=2 
	}
}

if(arbvar==FALSE){
	w=1
	while(w==1){
	if(is.null(x)){
	H0=normalmixEM(x=y,arbmean=TRUE,arbvar=FALSE,mu=mu, sigma=sigma, lambda=lambda,...)
	k=length(H0$lambda)
#	H1=normalmixEM(x=y,arbmean=TRUE,arbvar=TRUE,lambda=H0$lambda,mu=H0$mu,sigma=rep(H0$sigma,k),...)
	H1=normalmixEM(x=y,arbmean=TRUE,arbvar=TRUE,lambda=H0$lambda,mu=H0$mu,sigma=NULL,...)
	D=2*(H1$loglik-H0$loglik)
	df=k-1
	alpha=1-pchisq(D,df=df)
	} else{
	H0=regmixEM(y=y,x=x,arbmean=TRUE,arbvar=FALSE,beta=beta, sigma=sigma, lambda=lambda,...)
	k=length(H0$lambda)
#	H1=regmixEM(y=y,x=x,arbmean=TRUE,arbvar=TRUE,lambda=H0$lambda,beta=H0$beta,sigma=rep(H0$sigma,k),...)
	H1=regmixEM(y=y,x=x,arbmean=TRUE,arbvar=TRUE,lambda=H0$lambda,beta=H0$beta,sigma=NULL,...)
	D=2*(H1$loglik-H0$loglik)
	df=k-1
	alpha=1-pchisq(D,df=df)
	}
	if(D<0){
	w=1 
	mu=NULL
	sigma=NULL
	lambda=NULL
	} else w=2 
	}
}
a=list(chi.sq=D, df=df, p.value=alpha)
a
}
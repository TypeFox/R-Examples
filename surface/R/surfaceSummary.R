surfaceSummary <-function(fwd = NULL, bwd = NULL)
{

	if(is.null(fwd)) obj<-bwd
	if(is.null(bwd)) obj<-fwd

	if(!is.null(fwd)&!is.null(bwd)){
		obj<-c(fwd[-length(fwd)],bwd)
	}

	k<-length(obj)
	lnls<-sapply(obj,function(x){
			sapply(x$fit,function(y)summary(y)$loglik)
		})
	if(class(lnls)!="matrix")lnls<-matrix(lnls,nrow=1,dimnames=list(names(lnls)[1],NULL))
	aics<-sapply(obj,function(x)x$aic);names(aics)<-1:k
	n_regimes_seq<-sapply(obj,function(x)x$n_regimes)
	obj<-obj[[k]]
	shifts<-obj$savedshifts
	n_regimes<-obj$n_regimes
	obj<-obj$fit
	alpha<-sapply(obj,function(x)summary(x)$alpha)
	phylhalflife<-log(2)/alpha
	sigma_squared<-sapply(obj,function(x)summary(x)$sigma.squared)
	theta<-sapply(obj,function(x)summary(x)$optima[[1]])

	list(n_steps=k,lnls=lnls,n_regimes_seq=n_regimes_seq,aics=aics, shifts=shifts,n_regimes=n_regimes, alpha=alpha,phylhalflife=phylhalflife, sigma_squared=sigma_squared,theta=theta)
	}

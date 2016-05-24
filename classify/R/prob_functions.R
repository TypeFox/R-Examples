
rasch <- function(theta=NULL,delta=NULL){
	if(!is.numeric(theta)){
		stop("Person abilities need to be specified as a numeric vector.")		
	}
	if(!is.numeric(delta)){
		stop("Item difficulties need to be specified as a numeric vector.")		
	}
	p <- length(delta)
	probs <- sapply(seq(p), function(j) t(sapply(theta, function(th) exp(th-delta[j])/(1+exp(th-delta[j])))))
	return(probs)
}

tpl <- function(theta=NULL,beta=NULL,alpha=NULL){
	if(!is.numeric(theta)){
		stop("Person abilities need to be specified as a numeric vector.")		
	}
	if(!is.numeric(beta)){
		stop("Item difficulties need to be specified as a numeric vector.")		
	}
	if(!is.numeric(alpha)){
		stop("Item discrimination need to be specified as a numeric vector.")		
	}
	p <- length(beta)
	probs <- sapply(seq(p), function(j) t(sapply(theta, function(th) exp(alpha[j]*(th-beta[j]))/(1+exp(alpha[j]*(th-beta[j]))))))
	return(probs)
}

thpl <- function(theta=NULL,beta=NULL,alpha=NULL,eta=NULL){
	if(!is.numeric(theta)){
		stop("Person abilities need to be specified as a numeric vector.")		
	}
	if(!is.numeric(beta)){
		stop("Item difficulties need to be specified as a numeric vector.")		
	}
	if(!is.numeric(alpha)){
		stop("Item discrimination need to be specified as a numeric vector.")		
	}
	if(!is.numeric(eta)){
		stop("Item guessing parameters need to be specified as a numeric vector.")		
	}
	
	p <- length(beta)
	probs <- probs <- sapply(seq(p), function(j) t(sapply(theta, function(th) eta[j]+((1-eta[j]) * (exp(alpha[j]*(th-beta[j]))/(1+exp(alpha[j]*(th-beta[j]))))))))
	return(probs)
}

pcm <- function(theta=NULL,delta=NULL,n=NULL){
	if(!is.numeric(theta)){
		stop("Person abilities need to be specified as a numeric vector.")		
	}
	if(!is.numeric(delta)){
		stop("Item thresholds need to be specified as a numeric vector.")		
	}
	if(!is.numeric(n)){
		stop("Number of item categories must be specified")		
	}
	if(length(delta)!=n-1){
		stop("Delta thresholds should number one less than number of item categories")		
	}
	denor <- sum(exp(cumsum(theta-delta)))
	ps <- exp(cumsum(theta-delta))/(1+denor)
	z.p <- 1 - sum(ps)	
	ps <- c(z.p,ps,rep(0,n-length(ps)-1))
	return(ps)
}

gpcm <- function(theta=NULL,alpha=NULL,delta=NULL,n=NULL){
	if(!is.numeric(theta)){
		stop("Person abilities need to be specified as a numeric vector.")		
	}
	if(!is.numeric(delta)){
		stop("Item thresholds need to be specified as a numeric vector.")		
	}
	if(!is.numeric(alpha)){
		stop("Item discrimination needs to be specified as a numeric vector.")		
	}
	if(!is.numeric(n)){
		stop("Number of item categories must be specified")		
	}
	if(length(delta)!=n-1){
		stop("Delta thresholds should number one less than number of item categories")		
	}
	denor <- sum(exp(cumsum(alpha*(theta-delta))))
	ps <- exp(cumsum(alpha*(theta-delta)))/(1+denor)
	z.p <- 1 - sum(ps)	
	ps <- c(z.p,ps,rep(0,n-length(ps)-1))
	return(ps)
}



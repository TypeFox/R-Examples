# GIVEN GAMMA COMPUTE BETA, GAMMA, LOGLIK
# ALLOW FOR NULL Z OR NULL X

optim.lrt.lexpit <- function(beta.init,gamma.init,Y,X,Z,w,...){
	
	# NEED TO MATCH ARGUMENTS FOR OPTIM/NLM
	
LL <- function(Y,p,w){
		l <- w*(Y*logit(p)+log(1-p))
	-sum(l)
}

	if(is.null(Z)){ # BLM ONLY		
		beta.optim <- blm.optim(Y,X,w,beta.init)
		l <- LL(Y,X%*%beta.optim$par,w)
		gamma.optim <- numeric()	
	}
	else if(is.null(X)){ # LOGISTIC ONLY
		gamma.offset <- 0
		gamma.optim <- stage2.optim(Y,Z,gamma.offset,w,gamma.init)$est
		l <- LL(Y,expit(Z%*%gamma.optim),w)
		beta.optim <- list(par=numeric(), barrier.value=numeric())
	}
	else{
		
		beta.offset <- expit(Z%*%gamma.init)
		beta.optim <- stage1.optim(Y,X,w,beta.offset,beta.init)

		gamma.offset <- X%*%beta.optim$par
		gamma.optim <- stage2.optim(Y,Z,gamma.offset,w,gamma.init)$est
		
		l <- LL(Y,gamma.offset+expit(Z%*%gamma.optim),w)
	}
	
list(beta=beta.optim$par,gamma=gamma.optim,loglik=l,barrier.value=beta.optim$barrier.value)	
}


optim.lexpit <- function(beta.init,gamma.init,Y,X,Z,w,...){
	
	# NEED TO MATCH ARGUMENTS FOR OPTIM/NLM
	
	LL <- function(Y,p,w){
		l <- w*(Y*logit(p)+log(1-p))
	-sum(l)
}

	beta.offset <- expit(Z%*%gamma.init)
	beta.optim <- stage1.optim(Y,X,w,beta.offset,beta.init)
		
	gamma.offset <- X%*%beta.optim$par
	gamma.optim <- stage2.optim(Y,Z,gamma.offset,w,gamma.init)$est
		
	l <- LL(Y,gamma.offset+expit(Z%*%gamma.optim),w)
	
list(beta=beta.optim$par,gamma=gamma.optim,loglik=l,barrier.value=beta.optim$barrier.value)	
}

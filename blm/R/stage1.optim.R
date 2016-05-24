stage1.optim <- function(Y,X,w,offset,beta.init,...){

loglik <- function(beta,Y,X,w,offset){
	
	p <- X%*%beta+offset
	l <- w*(Y*logit(p)+log(1-p))

-sum(l)
}

gradient <- function(beta,Y,X,w,offset){
	
	p <- X%*%beta+offset
	C <- w*(Y*1/(p*(1-p))-1/(1-p))
	C <- matrix(C,nrow(X),ncol(X))

-colSums(C*X)
}

# DEFINE CONSTRAINTS
lexpit.constraints <- function(X,offset){
		
	# CONSTRAINTS OF FORM
	# X%*%beta-c >= 0
	list(
		U = rbind(X,-X),
		C = c(-offset,-(1-offset))	
	)
}

	constr <- lexpit.constraints(X,offset)

	results <- constrOptim(beta.init,f=loglik,grad=gradient,
						Y = Y, offset = offset, X = X, w = w,
					 	ui = constr$U,ci = constr$C,...)

results
}


blm.optim <- function(Y,X,w,beta.init,...){

loglik <- function(beta,Y,X,w){
	
	p <- X%*%beta
	l <- w*(Y*logit(p)+log(1-p))

-sum(l)
}

gradient <- function(beta,Y,X,w){
	
	p <- X%*%beta
	C <- w*(Y*1/(p*(1-p))-1/(1-p))
	C <- matrix(C,nrow(X),ncol(X))

-colSums(C*X)
}

# DEFINE CONSTRAINTS
lexpit.constraints <- function(X){
		
	# CONSTRAINTS OF FORM
	# X%*%beta-c >= 0
	list(
		U = rbind(X,-X),
		C = c(rep(0,nrow(X)),rep(-1,nrow(X)))	
	)
}

	constr <- lexpit.constraints(X)

	results <- constrOptim(beta.init,f=loglik,grad=gradient,
						Y = Y, X = X, w = w,
					 	ui = constr$U,ci = constr$C,...)

results
}

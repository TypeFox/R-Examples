vcov.blm.revised.bigmatrix <- function(beta,Y,X,w){
	
	x.params <- ncol(X)
	p <- X%*%beta
	A <- 1/nu(p)*(Y-p)
	nu.A <- (A+1)/(nu(p))
	tryDiag <- tryCatch(diag(as.numeric(w*nu.A)),error=function(x)NA)
	if(!is.matrix(tryDiag)){
		# COLUMN BY COLUMN PROCESSING
		W <- as.numeric(w*nu.A)
		XW <- X
		for(j in 1:x.params){
			XW[,j] <- X[,j]*W
		}
		HX <- -t(X)%*%XW
	}
	else{
		HX <- matrix(-t(X)%*%tryDiag%*%X,x.params,x.params)
	}
-HX
}

pearson.deviates <- function(object){
	
	p <- predict(object)
	w <- object@weights

sqrt(w)*(object@y-p)/sqrt(p*(1-p))
}

C <- function(object){
	
	X2 <- pearson.deviates(object)
	h <- leverage(object)

X2^2*h/(1-h)^2	
}

blm.leverage <- function(object){
	
	p <- predict(object)
	w <- object@weights
	X <- model.matrix(object@formula, object@data)
	V <- object@vcov
	
	nu <- p*(1-p)
	W <- diag(as.numeric(w/nu))
	HAT <- W%*%(X%*%V%*%t(X))

diag(HAT)
}

lexpit.leverage <- function(object){
	
	X <- cbind(model.matrix(object@formula.linear, object@data)[,-1])
	Z <- model.matrix(object@formula.expit, object@data)
	
	p.beta <- X%*%object@coef.linear
	p.gamma <- Z%*%object@coef.expit	
	w <- object@weights
	nu.beta <- p.beta*(1-p.beta)
	nu.gamma <- p.gamma*(1-p.gamma)
	
	W.BETA <- diag(as.numeric(w/nu.beta))
	W.GAMMA <- diag(as.numeric(w*nu.gamma))
	
	HAT.BETA <- W.BETA%*%(X%*%object@vcov.linear%*%t(X))
	HAT.GAMMA <- W.GAMMA%*%(Z%*%object@vcov.expit%*%t(Z))
	
diag(HAT.BETA+HAT.GAMMA)
}

vcov.blm.revised <- function(beta,Y,X,w){
	
	x.params <- ncol(X)
	p <- X%*%beta
	A <- 1/nu(p)*(Y-p)
	nu.A <- (A+1)/(nu(p))
	W <- diag(as.numeric(w*nu.A))
	HX <- matrix(-t(X)%*%W%*%X,x.params,x.params)
	
-HX
}


vcov.lexpit.revised <- function(beta,gamma,Y,X,Z,w){
	
	x.params <- ncol(X)
	z.params <- ncol(Z)
	n.params <- x.params+z.params
	
	p.beta <- X%*%beta
	p.gamma <- expit(Z%*%gamma)
	p <- p.beta+p.gamma
	A <- 1/nu(p)*(Y-p)
	nu.A <- (A+1)/(nu(p))
	
	W.beta <- diag(as.numeric(w*nu.A))
	W.gamma <- diag(as.numeric(w*ddexpit(p.gamma)*A-dexpit(p.gamma)^2*nu.A))
	W.beta.gamma <- diag(as.numeric(w*dexpit(p.gamma)*nu.A))
	
	HX <- -t(X)%*%W.beta%*%X
	HXZ <- -t(X)%*%W.beta.gamma%*%Z
	HZ <- t(Z)%*%W.gamma%*%Z
	
	H <- matrix(0,n.params,n.params)
	
	H[1:x.params,1:x.params] <- HX
	H[1:x.params,(x.params+1):(x.params+z.params)] <- HXZ
	H[(x.params+1):(x.params+z.params),1:x.params] <- t(HXZ)
	H[(x.params+1):(x.params+z.params),(x.params+1):(x.params+z.params)] <- HZ
	
-H
}

vcov.lexpit.revised.bigmatrix <- function(beta,gamma,Y,X,Z,w){
	
	x.params <- ncol(X)
	z.params <- ncol(Z)
	n.params <- x.params+z.params
	
	p.beta <- X%*%beta
	p.gamma <- expit(Z%*%gamma)
	p <- p.beta+p.gamma
	A <- 1/nu(p)*(Y-p)
	nu.A <- (A+1)/(nu(p))
	
	tryDiag <- tryCatch(diag(as.numeric(w*nu.A)),error=function(x)NA)
	if(!is.matrix(tryDiag)){
	 	W.beta <- as.numeric(w*nu.A)
	 	W.gamma <- as.numeric(w*ddexpit(p.gamma)*A-dexpit(p.gamma)^2*nu.A)
		W.beta.gamma <- as.numeric(w*dexpit(p.gamma)*nu.A)
		XW <- X
		ZW <- Z
		XZW <- Z
		
		for(j in 1:x.params){
			XW[,j] <- X[,j]*W.beta
		}
		
		for(j in 1:z.params){
			ZW[,j] <- Z[,j]*W.gamma
			XZW[,j] <- Z[,j]*W.beta.gamma
		}
	
		HX <- -t(X)%*%XW
		HXZ <- -t(X)%*%XZW
		HZ <- t(Z)%*%ZW
			
	}
	else{
		W.gamma <- diag(as.numeric(w*ddexpit(p.gamma)*A-dexpit(p.gamma)^2*nu.A))
		W.beta.gamma <- diag(as.numeric(w*dexpit(p.gamma)*nu.A))
	
		HX <- -t(X)%*%tryDiag%*%X
		HXZ <- -t(X)%*%W.beta.gamma%*%Z
		HZ <- t(Z)%*%W.gamma%*%Z
	}
	
	H <- matrix(0,n.params,n.params)
	
	H[1:x.params,1:x.params] <- HX
	H[1:x.params,(x.params+1):(x.params+z.params)] <- HXZ
	H[(x.params+1):(x.params+z.params),1:x.params] <- t(HXZ)
	H[(x.params+1):(x.params+z.params),(x.params+1):(x.params+z.params)] <- HZ
	
-H
}

vcov.blm.revised.strata <- function(beta,Y,X,w,strata){

	p <- ncol(X)
	n <- nrow(X)
	H <- matrix(0,p,p)
	
	for(i in levels(strata)){
		n.strata <- sum(strata==i)
		multiplier <- (n-1)/(n-p)*(n.strata/(n.strata-1))
		H <- H+multiplier*vcov.blm.revised.bigmatrix(beta,cbind(Y[strata==i,]),cbind(X[strata==i,]),cbind(w[strata==i,]))
	}

H
}

vcov.lexpit.revised.strata <- function(beta,gamma,Y,X,Z,w,strata){

	p <- ncol(X)+ncol(Z)
	n <- nrow(X)
	H <- matrix(0,p,p)
	
	for(i in levels(strata)){
		n.strata <- sum(strata==i)
		multiplier <- (n-1)/(n-p)*(n.strata/(n.strata-1))
		H <- H+multiplier*vcov.lexpit.revised(beta,gamma,
					cbind(Y[strata==i,]),cbind(X[strata==i,]),cbind(Z[strata==i,]),cbind(w[strata==i,]))
	}

H
}

# USE OF NLM PROCEDURE
stage2.optim <- function(Y,Z,offset,w,gamma.init,...){
	
stage2.f.optim <- function(gamma,Y,Z,offset,w){
	
loglik <- function(gamma,Y,Z,offset,w){
	p <- expit(Z%*%gamma)+offset
	l <- w*(Y*logit(p)+log(1-p))
-sum(l)
}

gradient <- function(gamma,Y,Z,offset,w){
	p <- expit(Z%*%gamma)+offset
	C <- w*dexpit(Z%*%gamma)*(Y*1/(p*(1-p))-1/(1-p))
	C <- matrix(C,nrow(Z),ncol(Z))
-colSums(C*Z)
}

hessian <- function(gamma,Y,Z,offset,w){
	
	p <- expit(Z%*%gamma)+offset
	C1 <- w*dexpit(Z%*%gamma)^2*(Y*1/(p*(1-p))^2*(1-2*p)+1/(1-p)^2)
	C2 <- w*ddexpit(Z%*%gamma)*(Y*1/(p*(1-p))-1/(1-p))
	ZZ <- apply(Z,1,function(z)outer(z,z))
	C <- matrix(C2-C1,ncol(Z)^2,nrow(Z),byrow=T)
	
-matrix(rowSums(C*ZZ),ncol(Z),ncol(Z))
}

	res <- loglik(gamma,Y,Z,offset,w)
	attr(res,"gradient") <- gradient(gamma,Y,Z,offset,w)
	
res
}

loglik <- function(gamma,Y,Z,offset,w){
	p <- expit(Z%*%gamma)+offset
	l <- w*(Y*logit(p)+log(1-p))
-sum(l)
}

fscale <- loglik(gamma.init,Y,Z,offset,w)

fit <- nlm(stage2.f.optim,gamma.init,Y=Y,Z=Z,offset=offset,w=w,print.level=0,
				check.analyticals=FALSE,typsize=gamma.init,
				fscale=fscale,...)

fit
}


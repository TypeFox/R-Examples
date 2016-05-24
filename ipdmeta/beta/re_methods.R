# TRIAL SPECIFIC DATA
Cpattern <- function(n,a,b){
	
	I <- diag(n)
	J <- rep(1,n)%*%t(rep(1,n))
	
a*I+b*J
}

solveCpattern <- function(n,a,b){
		Cpattern(n,1/a,-b/(a*(a+n*b)))
}

gfactor <- function(x,n,sigma2){x/(sigma2+n*x)}

W <- function(nt,nc,sigma2,D){
	
	# OBTAIN SUBJECT SPECIFIC INVERSE VARIANCE GIVEN SAMPLE SIZES
	# D MATRIX
	
	fe <- 1/sigma2*nc*(D[1,2]+D[1,1])^2*(1-nc*gfactor(D[1,1],nc,sigma2))
	fh <- 1/sigma2*nt*(D[1,2]+D[1,1])^2*(1-nt*gfactor(sum(D),nt,sigma2))
	E <- solveCpattern(nt,sigma2,sum(D)-fe)
	H <- solveCpattern(nc,sigma2,D[1,1]-fh)
	
	b1 <- b3 <- 1/sigma2
	b2 <- 1/sigma2*gfactor(sum(D),nt,sigma2)
	b4 <- 1/sigma2*gfactor(D[1,1]-fh,nc,sigma2)
	F <- -(D[1,1]+D[1,2])*(b1*b3-b1*b4*nc-b2*b3*nt+nc*nt*b2*b4)*(rep(1,nt)%*%t(rep(1,nc)))
	
rbind(cbind(E,F),cbind(t(F),H))
}

XW <- function(n,sigma2,D){
		
	# FOR BINARY CASE
	nt1 <- n[1]
	nt0 <- n[2]
	nc1 <- n[3]
	nc0 <- n[4]
	
	nt <- nt1+nt0
	nc <- nc1+nc0
	
	Wmat <- W(nt,nc,sigma2,D)
	
	a11 <- Wmat[1,2]
	a10 <- Wmat[1,1]-a11
	a2 <- Wmat[(nt+1),1]
	a31 <- Wmat[(nt+1),(nt+2)]
	a30 <- Wmat[(nt+1),(nt+1)]-a31
	
	X1 <- c(a10+nt*a11+a2*nc,a10+nt*a11+a2*nc,nt*a2+nc*a31+a30,nt*a2+nc*a31+a30)
	X2 <- c(a10+nt*a11,a10+nt*a11,nt*a2,nt*a2)
	X3 <- c(a10+nt1*a11+a2*nc1,nt1*a11+a2*nc1,nt1*a2+nc1*a31+a30,nt1*a2+nc1*a31)
	X4 <- c(a10+nt1*a11,nt1*a11,nt1*a2,nt1*a2)
	
rbind(X1,X2,X3,X4)
}


XWX <- function(n,sigma2,D){
	
	# FOR BINARY CASE
	nt1 <- n[1]
	nt0 <- n[2]
	nc1 <- n[3]
	nc0 <- n[4]
		
	nt <- nt1+nt0
	nc <- nc1+nc0
	
	Wmat <- W(nt,nc,sigma2,D)
	
	a11 <- Wmat[1,2]
	a10 <- Wmat[1,1]-a11
	a2 <- Wmat[(nt+1),1]
	a31 <- Wmat[(nt+1),(nt+2)]
	a30 <- Wmat[(nt+1),(nt+1)]-a31
	
	X1 <- c(a10+nt*a11+a2*nc,a10+nt*a11+a2*nc,nt*a2+nc*a31+a30,nt*a2+nc*a31+a30)
	X2 <- c(a10+nt*a11,a10+nt*a11,nt*a2,nt*a2)
	X3 <- c(a10+nt1*a11+a2*nc1,nt1*a11+a2*nc1,nt1*a2+nc1*a31+a30,nt1*a2+nc1*a31)
	X4 <- c(a10+nt1*a11,nt1*a11,nt1*a2,nt1*a2)
	
	XWX1 <- c(X1[1]*nt+X1[3]*nc,X1[1]*nt,X1[1]*nt1+X1[3]*nc1,X1[1]*nt1)
	XWX2 <- c(X1[1]*nt,X2[1]*nt,X2[1]*nt1+X2[3]*nc1,X2[1]*nt1)
	XWX3 <- c(XWX1[3],XWX2[3],X3[1]*nt1+X3[3]*nc1,X3[1]*nt1)
	XWX4 <- c(XWX1[4],XWX2[4],X3[1]*nt1,X4[1]*nt1)
	
rbind(XWX1,XWX2,XWX3,XWX4)
}


XWY <- function(y,n,sigma2,D){
		
	XWmat <- XW(n,sigma2,D)
	N <- diag(n)
	
XWmat%*%N%*%y
}

DZW <- function(nt,nc,sigma2,D){
	
	Wmat <- W(nt,nc,sigma2,D)
	
	a11 <- Wmat[1,2]
	a10 <- Wmat[1,1]-a11
	a2 <- Wmat[(nt+1),1]
	a31 <- Wmat[(nt+1),(nt+2)]
	a30 <- Wmat[(nt+1),(nt+1)]-a31
	
	DZW1 <- c((a10+a11*nt)*(D[1,1]+D[1,2])+D[1,1]*a2*nc,a2*nt*(D[1,1]+D[1,2])+D[1,1]*(a30+nc*a31))
	DZW2 <- c((a10+a11*nt)*(D[2,2]+D[1,2])+D[1,2]*a2*nc,a2*nt*(D[2,2]+D[1,2])+D[1,2]*(a30+nc*a31))
	
rbind(DZW1,DZW2)
}

ranef.meta <- function(y,n,sigma2,D,beta){
	
	mu <- c(sum(beta),sum(beta[1:2]),sum(beta[c(1,3)]),beta[1]) #ASSUME BINARY
	R <- y-mu
	N <- diag(n)
	D <- DZW(sum(n[1:2]),sum(n[3:4]),sigma2,D)  #ASSUME BINARY
	D <- cbind(D[,1],D[,1],D[,2],D[,2])

D%*%N%*%R
}

resid.meta <- function(s2,y,b,beta,n){
	
	Z <- matrix(c(1,1,1,1,1,1,0,0),4,2)
	X <- cbind(Z,Z)
	X[,3] <- c(1,0,1,0)
	X[,4] <- c(1,0,0,0)	
	mu <- X%*%beta+Z%*%b

sum((n-1)*s2+n*y^2-2*n*mu*y+n*mu^2)
}



lme.beta <- function(Y,N,sigma2,D){
	
	# DATA MATRIX OF TRIAL BY (TREATMENT-COVARIATE)
	
	XWXs <- sapply(1:nrow(Y),function(i){
		XWX(N[i,],sigma2,D)
	})
	
	solveXWX <- matrix(rowSums(XWXs),sqrt(nrow(XWXs)),sqrt(nrow(XWXs)))
	solveXWX <- solve(solveXWX)
	
	XYs <- sapply(1:nrow(Y),function(i){
		XWY(Y[i,],N[i,],sigma2,D)
	})
	
	XYs <- rowSums(XYs)

solveXWX%*%XYs
}

lm.meta <- function(Y,N,S2,initial=NULL,iter=10){
	
	if(is.null(initial)){
		initial <- list(sigma2=1,D=diag(2))
	}
	
lm.meta.loop <- function(Y,N,S2,Var){
	
	sigma2 <- Var[[1]]
	D <- Var[[2]]
	
	beta <- lme.beta(Y,N,sigma2,D)
	
	Y.list <- lapply(1:nrow(Y),function(i)Y[i,])
	N.list <- lapply(1:nrow(N),function(i)N[i,])
	S2.list <- lapply(1:nrow(S2),function(i)S2[i,])
	
	b <- mapply(ranef.meta,y=Y.list,n=N.list,
					MoreArgs=list(sigma2=sigma2,D=D,beta=beta),
					SIMPLIFY=FALSE)
	
	r <- mapply(resid.meta,s2=S2.list,y=Y.list,b=b,n=N.list,
						MoreArgs=list(beta=beta))

	sigma2 <- sum(r,na.rm=T)/sum(N)
	D <- matrix(rowSums(sapply(b,function(x)outer(x,x))),2,2)/length(b)

list(
	coef.fixed=beta,
	coef.random=b,
	cov.resid=sigma2,
	cov.random=D
)
}

	results <- list(lm.meta.loop(Y,N,S2,initial))
	
	for(i in 2:iter){
		results[[i]] <- lm.meta.loop(Y,N,S2,list(results[[i-1]]$cov.resid,
											results[[i-1]]$cov.random))
	}

	results <- list(
		coef.fixed = results[[iter]][[1]],
		coef.random = matrix(unlist(results[[iter]][[2]]),2),
		vcov.resid = results[[iter]][[3]],
		vcov.random = results[[iter]][[4]]
	)
	
	final.var <- mapply(XWX,n=N.list,
							MoreArgs=list(sigma2=results$vcov.resid,
										D=results$vcov.random))
										
	final.var <- matrix(rowSums(final.var),nrow(results$coef.fixed),nrow(results$coef.fixed))

	results$var <- solve(final.var)

results
}
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

Xblock <- function(n){

	X <- matrix(0,ncol=length(n)-1,nrow=sum(n))
	
	 for(i in 1:ncol(X)){ # (n-1) dummy variables
	   	start = 1
	   	if(i>1) start = start + sum(n[1:i-1])
	   	X[start:(start+n[i]-1),i] <- rep(1,n[i])
   }

X	
}

design.matrix <- function(n){ # CONSTRUCT WITHIN-STUDY DESIGN MATRIX
	
	intercept <- rep(1,sum(n))
	treatment <- rep(c(1,0),c(sum(n$trt),sum(n$ctrl)))
	covariates <- rbind(Xblock(n$trt),Xblock(n$ctrl))
	interactions <- rbind(Xblock(n$trt),Xblock(n$ctrl)*0)	

cbind(intercept,treatment,covariates,interactions)
}


XW <- function(a10,a11,a2,a30,a31,n){
	
	N <- function(n,constants){
		matrix(constants[-length(constants)],nrow=length(constants)-1,ncol=n)	
	}
	
	# n is a data.frame with trt and control, for each of K categories
	nt <- sum(n$trt)
	nc <- sum(n$ctrl)
	
	# TREATMENT, CONTROL
	intercept <- c(a10+nt*a11+nc*a2,nt*a2+nc*a31+a30)
	treatment <- c(a10+nt*a11,nt*a2)

		
   A1 <- rep(c(1,0),c(nt,nc))
   A2 <- A1
   
   A1 <- ifelse(A1==1,intercept[1],intercept[2])
   A2 <- ifelse(A2==1,treatment[1],treatment[2])
   
   XT <- t(Xblock(n$trt))
   XC <- t(Xblock(n$ctrl))
   XZT <- XT
   XZC <- XC
   
   T <- N(nt,n$trt)
   C <- N(nt,n$ctrl)
   XT <- a10*XT+a11*T+a2*C
   XZT <- a10*XZT+a11*T
  
  	T <- N(nc,n$trt)
	C <- N(nc,n$ctrl)
  	XC <- a30*XC+a2*T+a31*C
	XZC <- a2*T

  rbind(A1,A2,cbind(XT,XC),cbind(XZT,XZC))
}

quadratic.forms <- function(n,y.bar,sigma2,D){
	
	# CONSTRUCT STUDY-SPECIFIC XW, XWX MATRICES 
	# GIVEN SAMPLE SIZES AND VARIANCE COMPONENTS
	X <- design.matrix(n)
	W0 <- W(sum(n$trt),sum(n$ctrl),sigma2,D)
	XW0 <- t(X)%*%W0
	Y <- rep(c(y.bar$trt,y.bar$ctrl),c(n$trt,n$ctrl))
	
	list(
		XWY = XW0%*%Y, # CHECK FORM OF XW0 before processing
		XWX = XW0%*%X
		)
}

ipd_lme_ranef_study <- function(n,y.bar,fixef,sigma2,D){
	
	# STUDY INTERCEPT AND SLOPE
	
	X <- design.matrix(n)
	Z <- X[,1:2]
	W0 <- W(sum(n$trt),sum(n$ctrl),sigma2,D)
	DZW <- D%*%t(Z)%*%W0

	Y <- rep(c(y.bar$trt,y.bar$ctrl),c(n$trt,n$ctrl))
	R <- Y-X%*%fixef$fixef
	beta <-  DZW%*%R
	V <- D%*%t(Z)%*%(W0-W0%*%X%*%fixef$vcov%*%t(X)%*%W0)%*%Z%*%D
	
	list(
		ranef = beta,
		vcov.ranef = V
	)
}

ipd_lme_fixef <- function(n,y,sigma2,D){
	# LIST OF AGGREGATE DATA	
	qrs <- function(N,Y){
		quadratic.forms(N,Y,sigma2,D)	
	}
	
	components <- mapply(qrs,N=n,Y=y,SIMPLIFY=FALSE)
	
	XWYs <- sapply(components,function(x) x$XWY)
	XWXs <- sapply(components,function(x) x$XWX)
	
	XWY <- rowSums(XWYs)
	XWX <- rowSums(XWXs)
	p <- sqrt(length(XWX))
	XWX <- matrix(XWX,p,p)
	
	V <- solve(XWX)
	beta <- V%*%XWY

	labels <- ipd_lme_labels(n)
	row.names(beta) <- labels
	row.names(V) <- labels
	colnames(V) <- labels
	
list(
	fixef = beta,
	vcov.fixef = V
)

}

ipd_lme_ranef <- function(n,y,fixef,sigma2,D){
	mapply(ipd_lme_ranef_study,n=n,y.bar=y,
				  MoreArgs=list(fixef=fixef,sigma2=sigma2,D=D),SIMPLIFY=FALSE)
}

ipd_lme_labels <- function(N){
	
	labels <- c("Intercept","trt")
	x.labels <- row.names(N[[1]])[-nrow(N[[1]])]
	interaction.labels <- paste("trt",x.labels,sep=":")

c(labels,x.labels,interaction.labels)
}

ipd_lme_sigma2_study <- function(n,y,s2,beta,alpha){
	
	# STUDY INTERCEPT AND SLOPE
	N <- sum(c(n$trt,n$ctrl))
	y.bar <- sum(c(y$trt,y$ctrl)*c(n$trt,n$ctrl))/N	
	Y <- rep(y.bar,N)
	X <- design.matrix(n)
	Z <- X[,1:2]
	mu <- X%*%beta+Z%*%alpha
	SS <- t(Y-mu)%*%(Y-mu)
	Q <- (N-1)*s2-SS

Q
}

ipd_lme_sigma2 <- function(n,y,beta,alphas,s2){
	
	Qs <- mapply(ipd_lme_sigma2_study,n=n,y=y,alpha=alphas,s2=s2,MoreArgs=list(beta=beta))
	N <- sapply(n,function(x)sum(c(x$trt,x$ctrl)))
	N <- sum(N)
	
sum(Qs)/(N-length(beta))
}

ipd_lme_D <- function(alphas){
	
	Q <- sapply(alphas,function(x)outer(x,x))
	Q <- matrix(rowSums(Q),2,2)
	
Q/length(alphas)
}

initialize_vcov <- function(n, y, s2){
	
	# FOR RESIDUAL TAKE MEAN OF S2
	NT <- sapply(n,function(x)sum(x$trt))
	NC <- sapply(n,function(x)sum(x$ctrl))
	
	YT <- mapply(function(n,y){sum(n$trt*y$trt)/sum(n$trt)},n=n,y=y)
	YC <- mapply(function(n,y){sum(n$ctrl*y$ctrl)/sum(n$ctrl)},n=n,y=y)
	
	sigma2 <- sum((NT+NC-1)*s2)/(sum(NT+NC)-1)
	y.bar.ctrl <-  sum((NC-1)*(YC-mean(YC))^2)/sum(NC)
	y.bar.trt <- sum((NT-1)*(((YT-mean(YT))-(YC-mean(YC))))^2)/sum(NT)
	
	list(
		sigma2 = sigma2,
		D = diag(c(y.bar.ctrl,y.bar.trt))
		)	
}

ipd_lme_meta <- function(n, y, s2,max.iter=50,tol=1e-7){
	
	vcov.init <- initialize_vcov(n, y, s2)
	beta.trace <- NULL
	alpha.trace <- NULL
	sigma.trace <- NULL
	D.trace <- NULL
	convergence.sigma <- tol
	convergence.D <- tol
			
	# UPDATE FIXED POP AND STUDY EFFECTS
	for(j in 1:max.iter){

		beta <- ipd_lme_fixef(n,y,sigma2=vcov.init$sigma2,D=vcov.init$D)
		alpha <- ipd_lme_ranef(n,y,beta,sigma2=vcov.init$sigma2,D=vcov.init$D)
	
    	# UPDATE VARIANCE COMPONENTS
		alphas <- lapply(alpha,function(x)x$ranef)
		vcov.init$sigma2 <- ipd_lme_sigma2(n,y,s2,alpha=alphas,beta=beta$fixef)
		vcov.init$D <- ipd_lme_D(alphas)

		beta.trace <- cbind(beta.trace,beta$fixef)
		alpha.trace <- cbind(alpha.trace,unlist(alphas))
		sigma.trace <- cbind(sigma.trace,vcov.init$sigma2)
		D.trace <- cbind(D.trace,as.numeric(vcov.init$D))	
		
		
		if(j>1){
			convergence.sigma <- (sigma.trace[,ncol(sigma.trace)-1]-sigma.trace[,ncol(sigma.trace)])^2/apply(sigma.trace,1,median)
			convergence.D <- max((D.trace[,ncol(D.trace)-1]-D.trace[,ncol(D.trace)])^2/abs(apply(D.trace,1,median)	))	
		}
		
		if(convergence.sigma<tol&convergence.D<tol) break
		
	}
	
	list(
		   beta=beta,
		   alpha=alpha,
		   vcov = vcov.init,
		   beta.trace = beta.trace,
		   alpha.trace = alpha.trace,
		   sigma.trace = sigma.trace,
		   D.trace = D.trace,
		   n.iter = j,
		   max.iter = max.iter,
		   tol = tol
		   )
}


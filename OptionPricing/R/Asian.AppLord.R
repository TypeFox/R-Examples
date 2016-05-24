BS_A <- function(K,beta,gamma){
	if(K>0){
		m <- (beta-log(K))/gamma
		exp(beta+gamma^2/2)*pnorm(m+gamma)-K*pnorm(m)
	}
	else
	exp(beta+gamma^2/2)
}	


cond2M_A <- function(z,T=1,d=12,r=0.05,sigma=0.1,S0=100,C,evalC=TRUE,evalm=TRUE){
	# mean and variance of A conditional on G

	dt <- T/d
	if(evalC==TRUE) C <- Covmat_A(T,d,sigma)

	logS0 <- log(S0)

	# Mean and sd of logGA
	mus <- logS0+(r-sigma^2/2)*dt*(d+1)/2
	sigmas <- sigma/d*sqrt(dt*d*(d+1)*(2*d+1)/6)

	# Expectation and variance of logSti
	muls <- function(i) logS0 + (r-sigma^2/2)*dt*i
	sigma2ls <- function(i) sigma^2*dt*i

	# Covariance of log(GA) and logSti 
	covi <- function(i) sigma^2*dt/d * (i*(i+1)/2+(d-i)*i)
 
	muhat <- function(i) muls(i) + covi(i)*z/sigmas

	if(evalm==TRUE) meanA <- (1/d)* sum(exp(muhat(1:d)+.5*diag(C)))

	muhatvec <- muhat(1:d)
	diagC <- diag(C)
	vec <- muhatvec+.5*diagC
	evec <- exp(vec)
	Sum <- sum(evec)
	mat <- matrix(vec,d,d,byrow=T)
	varA <- sum(evec*(rowSums(exp(mat+C))-Sum))
	varA <- varA/d^2
	if(evalm==TRUE) c(meanA,varA) else varA
}

Covmat_A <- function(T=1,d=12,sigma=0.1){
	dt <- T/d
	# Covariance matrix of logSti
	C <- matrix(0,d,d)
	for(i in 1:d) C[i,i:d] <- sigma^2*dt*i
	for(i in 2:d) C[i,1:i] <- C[1:i,i]
	
	# Covariance with logG
	w <- rep(1,d)/d
	Cw <- C %*% w
	#print(Cw)
	C - (Cw %*% t(Cw))/(w%*%C%*%w)[1,1]
}


evalECV_A <- function(T=1,d=12,K=100,r=0.05,sigma=0.1,S0=100){
	# Closed form solution of the expectation of the new CV
	# taken from Curran(1994)
	# see formula (8) in Dingec and Hormann (2013)

	dt <- T/d
	logS0 <- log(S0)
	logK <- log(K)

	# Mean and sd of logGA
	mus <- logS0+(r-sigma^2/2)*dt*(d+1)/2
	sigmas <- sigma/d*sqrt(dt*d*(d+1)*(2*d+1)/6)
	k <- (logK-mus)/sigmas
	I2 <- K*pnorm(-k)	

	# Expectation and variance of logSti
	muls <- function(i) logS0 + (r-sigma^2/2)*dt*i
	sigma2ls <- function(i) sigma^2*dt*i
	
	# Covariance of log(GA) and logSti 
	covi <- function(i) sigma^2*dt/d * (i*(i+1)/2+(d-i)*i)

	I1 <- (1/d)*sum(exp(muls(1:d)+sigma2ls(1:d)/2)*pnorm(-k+covi(1:d)/sigmas))
	
	exp(-r*T)*(I1-I2)
}
 
AsianCall_AppLord <- function(T=1,d=12,K=100,r=0.05,sigma=0.1,S0=100,all=TRUE){
	# Approximation of Lord (2006)
	# all ... if TRUE, approximation is given for the whole price
	dt <- T/d
	mus <- log(S0)+(r-sigma^2/2)*dt*(d+1)/2
	sigmas <- sigma/d*sqrt(dt*d*(d+1)*(2*d+1)/6)
	k <- (log(K)- mus)/sigmas
	C <- Covmat_A(T,d,sigma)
	q <- function(z){
		momvec <- cond2M_A(z,T,d,r,sigma,S0,C,evalC=F,evalm=TRUE)
		me <- momvec[1]- exp(mus+sigmas*z); 
		ve <- momvec[2]
		gamma <- sqrt(log(ve/me^2+1)) 
		beta <- log(me)-gamma^2/2
		#print(c(beta,gamma))	
		BS_A(K=K-exp(mus+sigmas*z),beta,gamma)*dnorm(z)
	}
	res <- exp(-r*T)*integrate(Vectorize(q),-Inf,k)[[1]]
	if(all) res <- res+evalECV_A(T,d,K,r,sigma,S0)
	res
}







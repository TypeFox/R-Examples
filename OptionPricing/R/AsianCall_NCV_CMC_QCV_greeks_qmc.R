
#library(EOP)

# CV+CMC method for greeks


# Expectation of the arithmetic average

#expectA <- function(S0=100,T=1,d=12,r=0.05) S0/d*sum(exp(r*(1:d)*T/d))
#expectA()
# strike price as a function of moneyness
#(1+0.0)*expectA(T=5,d=5)


# Expectation of main CV

evalECV <- function(T=1,d=12,K=100,r=0.05,sigma=0.1,S0=100,greeks=FALSE){
	# Closed form solution of the expectation of the new CV
	# taken from Curran(1994)
	# see formula (8) in Dingec and Hormann (2013)
	dt <- T/d
	varX <- d*(d+1)*(2*d+1)/6
	mus <- log(S0)+(r-sigma^2/2)*dt*(d+1)/2
	sigmas <- sigma/d*sqrt(dt*varX)
	k <- (log(K)-mus)/sigmas
	v <- (d:1)/sqrt(varX)
	a <- sigma*sqrt(dt)*cumsum(v)
	sum1 <- sum(exp(r*(1:d)*dt)*pnorm(-k+a))
	price <- (S0/d)*sum1-K*pnorm(-k)
	if(greeks){
		hk <- (S0/d)*sum(exp(a*k+r*(1:d)*dt-a^2/2))
		hdk <- (S0/d)*sum(a*exp(a*k+r*(1:d)*dt-a^2/2))
		delta <- (1/d)*sum1+(hk-K)*dnorm(k)/(S0*sigmas)
		#gamma <- (hdk+(hk-K)*(sigmas-k))*dnorm(k)/(S0*sigmas)^2
		gamma <- (2*hk*sigmas-hdk-(hk-K)*(sigmas-k))*dnorm(k)/(S0*sigmas)^2
		res <- exp(-r*T)*c(price,delta,gamma)
		names(res) <- c("price","delta","gamma")
	}
	else res <- exp(-r*T)*price
	res 
}
#evalECV()


#evalECV(greeks=TRUE)


# checking the formula with finite difference
#evalECV(T=5,d=5,sigma=0.5,K=(1-0.0)*expectA(T=5,d=5),S0=100,greeks=TRUE)
#h<- 1e-3
#y0 <- evalECV(T=5,d=5,sigma=0.5,K=(1-0.0)*expectA(T=5,d=5),S0=100)
#y1 <- evalECV(T=5,d=5,sigma=0.5,K=(1-0.0)*expectA(T=5,d=5),S0=100+h)
#y2 <- evalECV(T=5,d=5,sigma=0.5,K=(1-0.0)*expectA(T=5,d=5),S0=100-h)
## delta
#(y1-y2)/(2*h)
## gamma
#(y1-2*y0+y2)/h^2
#
#(y1-y2)/(2*h)-evalECV(T=5,d=5,sigma=0.5,K=(1-0.0)*expectA(T=5,d=5),S0=100,greeks=TRUE)[2]
#(y1-2*y0+y2)/h^2-evalECV(T=5,d=5,sigma=0.5,K=(1-0.0)*expectA(T=5,d=5),S0=100,greeks=TRUE)[3]



newtons_root <- function(x0,func,maxiter=10,tol=1e-15,rhs=TRUE){
	# Newton's method for root finding
	# x0...initial point
	# func ... function
	if(rhs==TRUE) sign <- -1 else sign <- 1
	mat <- func(x0)
	for (i in 1:maxiter){
		x1 <- x0 - mat[,1]/mat[,2]
		mat <- func(x1)
		error <- x1-x0
		if(max(abs(error))<tol) break;
		#if(max(abs(mat[,1]))<tol) break;
		#print(max(abs(mat[,1])))
		x0 <- x1
	}
	#print(max(abs(mat[,1])))
	cbind(x1,error,i)
}




findrootapp <- function(a,s,k,K,d,n){
	# First aorder apprximations for the root
	# Approximate root
	# Derivatives
	mat <- exp(a*k)*s
	fk <- colSums(mat)
	fdk <- colSums(a*mat)
	#f3dk <- colSums(a^3*exp(a*k)*s)

	# f is convex: as sum of convex functions
	# all derivatives are positive

	#print((1/d)*cbind(fdk,f2dk,f3dk))
	# 1st order app
	bapp1 <- k+(d*K-fk)/fdk
	alpha <- fk/d	
	beta <- fdk/d
	res <- k+ alpha/beta*log(K/alpha)
	res
}

findbcv <- function(K=100,T=1,d=12,r=0.05,sigma=0.1,S0=100){
	# Finds the lower bound of the integral
	
	dt <- T/d
	logS0 <- log(S0)

	# Mean and sd of logGA
	mus <- logS0+(r-sigma^2/2)*dt*(d+1)/2
	sigmas <- sigma/d*sqrt(dt*d*(d+1)*(2*d+1)/6)
	
	#Log(GA)
	#print(paste("GA=",exp(x)))

	# Expectation and variance of logSti
	muls <- function(i) logS0 + (r-sigma^2/2)*dt*i
	sigma2ls <- function(i) sigma^2*dt*i
	
	# Covariance of log(GA) and logSti 
	covi <- function(i) sigma^2*dt/d * (i*(i+1)/2+(d-i)*i)

	f <- function(z){
		x <- mus+sigmas*z
		(1/d)*sum(exp(muls(1:d)+covi(1:d)/sigmas^2*(x-mus)+(sigma2ls(1:d)-covi(1:d)^2/sigmas^2)/2))-K
	}
	uniroot(f,c(-10,10),tol=1e-12)$root
} 




evalLB <- function(K=100,T=1,d=12,r=0.05,sigma=0.1,S0=100,greeks=FALSE,all=FALSE){
	# Closed form solution for the lower bound of Curran (1994)
	# see formula (26) in Dingec and Hormann (2013) p. 430.
	# all ... if TRUE, lower bound is given for the whole price
	
	dt <- T/d
	varX <- d*(d+1)*(2*d+1)/6
	mus <- log(S0)+(r-sigma^2/2)*dt*(d+1)/2
	sigmas <- sigma/d*sqrt(dt*varX)
	k <- (log(K)-mus)/sigmas
	v <- (d:1)/sqrt(varX)
	a <- sigma*sqrt(dt)*cumsum(v)
	bcv <- findbcv(K,T,d,r,sigma,S0)
	sum1 <- sum(exp(r*(1:d)*dt)*(pnorm(k-a)-pnorm(bcv-a)))
	price <- (S0/d)*sum1-K*(pnorm(k)-pnorm(bcv))
	if(greeks){
		hk <- (S0/d)*sum(exp(a*k+r*(1:d)*dt-a^2/2))
		hdk <- (S0/d)*sum(a*exp(a*k+r*(1:d)*dt-a^2/2))
		k0 <- -1/(S0*sigmas)
		h0bcv <- (1/d)*sum(exp(a*bcv+r*(1:d)*dt-a^2/2))
		hpbcv <- (S0/d)*sum(a*exp(a*bcv+r*(1:d)*dt-a^2/2))
		bcv0 <- -h0bcv/hpbcv 
		delta <- (1/d)*sum1 - (hk-K)*dnorm(k)/(S0*sigmas)
		# gamma
		C1 <- 1/S0*(hk*dnorm(k)*k0-K*dnorm(bcv)*bcv0)
		C2 <- dnorm(k)/(S0*sigmas)^2*(hdk-hk*sigmas+(hk-K)*(sigmas-k))
		gamma <- C1+C2
		res <- exp(-r*T)*c(price,delta,gamma)
		names(res) <- c("price","delta","gamma")
	}
	else res <- exp(-r*T)*price
	if(all) res <- res + evalECV(T,d,K,r,sigma,S0,greeks) 
	res 
}
#evalLB()

#evalLB(greeks=TRUE)


# checking the greek formulas with finite difference

#evalLB(T=5,d=5,sigma=0.5,K=(1-0.0)*expectA(T=5,d=5),S0=100,greeks=TRUE)
#
#h<- 1e-2
#y0 <- evalLB(T=5,d=5,sigma=0.5,K=(1-0.0)*expectA(T=5,d=5),S0=100)
#y1 <- evalLB(T=5,d=5,sigma=0.5,K=(1-0.0)*expectA(T=5,d=5),S0=100+h)
#y2 <- evalLB(T=5,d=5,sigma=0.5,K=(1-0.0)*expectA(T=5,d=5),S0=100-h)
## delta
##(y1-y2)/(2*h)
## gamma
##(y1-2*y0+y2)/h^2
#
#(y1-y2)/(2*h)-evalLB(T=5,d=5,sigma=0.5,K=(1-0.0)*expectA(T=5,d=5),S0=100,greeks=TRUE)[2]
#(y1-2*y0+y2)/h^2-evalLB(T=5,d=5,sigma=0.5,K=(1-0.0)*expectA(T=5,d=5),S0=100,greeks=TRUE)[3]




eval_equad <- function(T=1,d=12,K=100,r=0.05,sigma=0.1,S0=100){
	# Expectation of the quadratic CV
	# see formula (27) in Dingec and Hormann (2013) p. 430.
	
	dt <- T/d
	sumSt <- sumlogSt <-0
	#logSt <- log(S0)

	varX <- d*(d+1)*(2*d+1)/6
	mus <- log(S0)+(r-sigma^2/2)*dt*(d+1)/2
	sigmas <- sigma/d*sqrt(dt*varX)
	k <- (log(K)-mus)/sigmas
	v <- (d:1)/sqrt(varX)
	a <- sigma*sqrt(dt)*cumsum(v)
	bcv <- findbcv(K,T,d,r,sigma,S0)
	mcv <- exp(a^2/2)*(pnorm(k-a)-pnorm(bcv-a))

	#(1/d)*colSums(mcv*s)-K*(pnorm(k)-pnorm(bcv))

	#term1 <- E[((1/d)*colSums(mcv*s))^2]
	#term2 <- -2*(1/d)*K*(pnorm(k)-pnorm(bcv)) * E[colSums(mcv*s)]
	term3 <- (K*(pnorm(k)-pnorm(bcv)))^2

	mu <- (r-sigma^2/2)*(1:d)*dt
	
	I <- diag(1, d, d)
	varcov <- I-(v%o%v)
	L <- lower.tri(I,TRUE)*1
	varcov <- L %*% varcov %*% t(L)

	vr <- diag(varcov)
	matvr <- matrix(vr,d,d)
	sigmas2 <- sigma^2*dt*(matvr + t(matvr) + 2* varcov)

	matmus <- matrix(mu,d,d)
	mus <- matmus + t(matmus)

	term2 <- sum(mcv*S0*exp(mu+(sigma^2*dt*vr)/2))
	term2 <- term2*(-2)*(1/d)*K*(pnorm(k)-pnorm(bcv))
	
	mmat <- matrix(mcv,d,d)
	term1 <- (1/d)^2*S0^2*sum(mmat*t(mmat)*exp(mus+sigmas2/2))

	#print(c(term1,term2,term3))
	term1+term2+term3
}

AsianCall_CMC_CV <- function(n=10^4,T=1,d=12,K=100,r=0.05,sigma=0.1,S0=100,
				maxiter=100,tol=1e-14,np=100,plotyn=FALSE){
	# New simulation algorithm for Asian call options 
	# new CV, conditional Monte Carlo and quadratic CVs
	# Algorithm 7 in Dingec and Hormann (2013) p. 430.
	# n ... number of repetitions for the simulation
	# tol ... tolerance of Newton's method
	# maxiter ...  max number of iterations for Newton's method
	# plotyn ... if TRUE, the scatter plot of y against psi is plotted.
	# returncs ... if TRUE, optimal regression coefficients are returned.
	# pilotrun .. if TRUE, regression coefficients are estimated in a pilot simulation run, otherwise, all sample is used.
	# np ... sample size of pilot simulation

	n <- n+np
	dt <- T/d
	sumSt <- sumlogSt <-0
	logSt <- log(S0)

	varX <- d*(d+1)*(2*d+1)/6
	mus <- log(S0)+(r-sigma^2/2)*dt*(d+1)/2
	sigmas <- sigma/d*sqrt(dt*varX)
	k <- (log(K)-mus)/sigmas

	# Simulate marginals conditional on Z
	v <- (d:1)/sqrt(varX)
	Z <- matrix(rnorm(n*d),d,n)
	zmatr <- Z - v%o%(v%*%Z)[1,]
	
	a <- sigma*sqrt(dt)*cumsum(v)
	s <- matrix(0,d,n)
	s[1,] <- S0*exp((r-sigma^2/2)*dt+sigma*sqrt(dt)*zmatr[1,])
	for(i in 2:d)
		s[i,] <- s[i-1,]*exp((r-sigma^2/2)*dt+sigma*sqrt(dt)*zmatr[i,])	
	
	# Newton-Raphson method to find the root: b(Z)

	# Function as an argument to Newtons method
	# Returns the function value and its derivative at x	
	fdf <- function(x) {
		mat <- exp(a%o%x)*s
		cbind((1/d)*colSums(mat)-K,(1/d)*colSums(a*mat))
	}
	# Approximate root as a starting point (1st order app)
	x0 <- findrootapp(a,s,k,K,d,n)
	# Newtons method
	nwroot <- newtons_root(x0,fdf,maxiter,tol)
	b <- nwroot[,1]

	# Computing conditional expectation
	m <- matrix(0,d,n)
	for(i in 1:d) m[i,] <- exp(a[i]^2/2)*(pnorm(k-a[i])-pnorm(b-a[i]))	
	y <- exp(-r*T)*((1/d)*colSums(m*s)-K*(pnorm(k)-pnorm(b)))
	
	# CV
	bcv <- findbcv(K,T,d,r,sigma,S0)
	mcv <- exp(a^2/2)*(pnorm(k-a)-pnorm(bcv-a))
	cecv <- (1/d)*colSums(mcv*s)-K*(pnorm(k)-pnorm(bcv))
	x1 <- exp(-r*T)* cecv
	x2 <- cecv^2

	if(plotyn==TRUE){
		plot(x1,y,cex=.5,xlab=expression(psi),ylab="Y")
		#lines(sort(cecv),sort(reg$fitted.values))
	}
	
	# Linear Regression
	cs <- lm(y[1:np]~x1[1:np]+x2[1:np])$coeff[-1]
	ex1 <- evalLB(K,T,d,r,sigma,S0,FALSE)
	ex2 <- eval_equad(T,d,K,r,sigma,S0)
	#print(var(y)/var(y-cs[1]*(x1-ex1)-cs[2]*(x2-ex2)))

	ycv <- y-cs[1]*(x1-ex1)-cs[2]*(x2-ex2)
	ycv <- ycv[-(1:np)] + evalECV(T,d,K,r,sigma,S0) 
	n <- n-np
	error <- 1.96*sd(ycv)/sqrt(n)
	res <- c(mean(ycv),error)
	names(res)<-c("result","error estimate")
	res
}

#AsianCall_CMC_CV()

#AsianCall_CVimpCMCCV(n = 10^4, T = 1, d = 12, K = 100, r = 0.05, sigma = 0.1, S0 = 100)


# Naive simulation of delta and gamma

AsianCall_naive_greeks <- function(n=10^4,T=1,d=12,K=100,r=0.05,sigma=0.2,S0=100){
	# Naive Simulation Algorithm for Asian option with arithmetic average
	# Algorithm 2 in Dingec and Hormann (2013), p.423
	# n ... number of repetitions for the simulation
	dt <- T/d
	sumSt <- 0
	logSt <- log(S0)

	for(i in 1:d){
		Z <- rnorm(n)
		logSt <- logSt+(r-sigma^2/2)*dt+sigma*sqrt(dt)*Z
		sumSt <- sumSt + exp(logSt)
		if(i==1){
			Z1 <- Z
			L1 <- sigma*sqrt(dt)
		}
	}
	A <- sumSt/d
	price <- exp(-r*T)*pmax(A-K,0)	
	delta <- exp(-r*T)*(A>K)*A/S0
	# mixed estimator (LR-PW)
	gamma <- exp(-r*T)*(A>K)*K*Z1/(S0^2*L1)
	
	error_price <- 1.96*sd(price)/sqrt(n)
	error_delta <- 1.96*sd(delta)/sqrt(n)
	error_gamma <- 1.96*sd(gamma)/sqrt(n)
	
	res <- rbind(c(mean(price),error_price),c(mean(delta),error_delta),c(mean(gamma),error_gamma))
	rownames(res)<-c("price","delta","gamma")
	colnames(res)<-c("result","error estimate")
	res
}

#AsianCall_naive_greeks()


# NCV+LR


AsianCall_NCV_LR_greeks <- function(n=10^4,T=1,d=12,K=100,r=0.05,sigma=0.1,S0=100){
	# Naive Simulation Algorithm for Asian option with arithmetic average
	# Algorithm 2 in Dingec and Hormann (2013), p.423
	# n ... number of repetitions for the simulation
	dt <- T/d
	sumSt <- sumlogSt <- 0
	logSt <- log(S0)

	for(i in 1:d){
		Z <- rnorm(n)
		logSt <- logSt+(r-sigma^2/2)*dt+sigma*sqrt(dt)*Z
		sumSt <- sumSt + exp(logSt)
		sumlogSt <- sumlogSt + logSt
		if(i==1){
			Z1 <- Z
			L1 <- sigma*sqrt(dt)
		}
	}
	A <- sumSt/d
	G <- exp(sumlogSt/d)
	
	y <- exp(-r*T)*pmax(A-K,0)*(G<K)
	ew <- evalECV(T,d,K,r,sigma,S0,TRUE)
	price <- y + ew[1] 
	delta <- y*Z1/(S0*L1) + ew[2]
	gamma <- y*(Z1^2-Z1*L1-1)/(S0*L1)^2 + ew[3]	
	
	error_price <- 1.96*sd(price)/sqrt(n)
	error_delta <- 1.96*sd(delta)/sqrt(n)
	error_gamma <- 1.96*sd(gamma)/sqrt(n)
	
	res <- rbind(c(mean(price),error_price),c(mean(delta),error_delta),c(mean(gamma),error_gamma))
	rownames(res)<-c("price","delta","gamma")
	colnames(res)<-c("result","error estimate")
	res
}

#AsianCall_NCV_LR_greeks()


# CMC + CV

# Expectation formulas for QCVs

evalEQCV <- function(T=1,d=12,K=100,r=0.05,sigma=0.1,S0=100,greeks=FALSE){
	# Expectation of the quadratic CV
	# see formula (27) in Dingec and Hormann (2013) p. 430.
	
	dt <- T/d
	#logSt <- log(S0)

	varX <- d*(d+1)*(2*d+1)/6
	mus <- log(S0)+(r-sigma^2/2)*dt*(d+1)/2
	sigmas <- sigma/d*sqrt(dt*varX)
	k <- (log(K)-mus)/sigmas
	v <- (d:1)/sqrt(varX)
	a <- sigma*sqrt(dt)*cumsum(v)
	bcv <- findbcv(K,T,d,r,sigma,S0)
	
	gamma <- (1/d)*exp(a^2/2)*(pnorm(k-a)-pnorm(bcv-a))
	eta <- pnorm(k)-pnorm(bcv)
	term3 <- (K*eta)^2

	mu <- (r-sigma^2/2)*(1:d)*dt
	I <- diag(1, d, d)
	varcov <- I-(v%o%v)
	L <- lower.tri(I,TRUE)*1
	varcov <- L %*% varcov %*% t(L)

	vr <- diag(varcov)
	H <- exp(mu+(sigma^2*dt*vr)/2)
	term2 <- -2*S0*sum(gamma*H)*K*eta
	
	matvr <- matrix(vr,d,d)
	sigmas2 <- sigma^2*dt*(matvr + t(matvr) + 2* varcov)
	matmus <- matrix(mu,d,d)
	mus <- matmus + t(matmus)
	M <- exp(mus+sigmas2/2)
	gmat <- matrix(gamma,d,d)
	term1 <- S0^2*sum(gmat*t(gmat)*M)

	#print(c(term1,term2,term3))
	res <- price <- term1+term2+term3
	
	if(greeks){
		k0 <- -1/(S0*sigmas)
		k00 <- 1/(S0^2*sigmas)
		h0bcv <- (1/d)*sum(exp(a*bcv+r*(1:d)*dt-a^2/2))
		hpbcv <- (S0/d)*sum(a*exp(a*bcv+r*(1:d)*dt-a^2/2))
		hppbcv <- (S0/d)*sum(a^2*exp(a*bcv+r*(1:d)*dt-a^2/2))
		bcv0 <- -h0bcv/hpbcv
		bcv00 <- -bcv0*(2/S0+ hppbcv* bcv0/hpbcv)
		#print(bcv00)
		gamma0 <- (1/d)*exp(a^2/2)*(dnorm(k-a)*k0-dnorm(bcv-a)*bcv0)
		gamma00 <-(1/d)*exp(a^2/2)*(dnorm(k-a)*(k00-k0^2*(k-a))-dnorm(bcv-a)*(bcv00-bcv0^2*(bcv-a)))
		eta0 <- dnorm(k)*k0-dnorm(bcv)*bcv0
		eta00 <- dnorm(k)*(k00-k0^2*k)-dnorm(bcv)*(bcv00-bcv0^2*bcv) 
		gmat <- matrix(gamma,d,d)
		gmat0 <- matrix(gamma0,d,d)
		gmat00 <- matrix(gamma00,d,d)

		# Delta
		mat1 <- 2*S0*gmat*t(gmat)+ S0^2*(gmat0*t(gmat)+gmat*t(gmat0))
		vec1 <- gamma*eta+S0*(gamma0*eta+gamma*eta0)
		Delta <- sum(mat1*M) - 2*K*sum(vec1*H) + 2*K^2*eta*eta0
		# Gamma
		mat2 <- 2*gmat*t(gmat) + 4*S0*(gmat0*t(gmat)+gmat*t(gmat0)) +S0^2*(gmat00*t(gmat)+2*gmat0*t(gmat0)+gmat*t(gmat00))
		vec2 <- 2*(gamma0*eta+gamma*eta0)+S0*(gamma00*eta+2*gamma0*eta0+gamma*eta00)
		Gamma <- sum(mat2*M) -2*K*sum(vec2*H) + 2*K^2*(eta0^2+eta*eta00)
		res <- c(price,Delta,Gamma)
		names(res) <- c("price","delta","gamma")	
	}
	res
}

#evalEQCV()
#evalEQCV(greeks=TRUE)


# checking the greek formulas with finite difference

#evalEQCV(T=5,d=5,sigma=0.5,K=(1-0.0)*expectA(T=5,d=5),S0=100,greeks=TRUE)
#
#h<- 1e-2
#y0 <- evalEQCV(T=5,d=5,sigma=0.5,K=(1-0.0)*expectA(T=5,d=5),S0=100)
#y1 <- evalEQCV(T=5,d=5,sigma=0.5,K=(1-0.0)*expectA(T=5,d=5),S0=100+h)
#y2 <- evalEQCV(T=5,d=5,sigma=0.5,K=(1-0.0)*expectA(T=5,d=5),S0=100-h)
## delta
##(y1-y2)/(2*h)
## gamma
##(y1-2*y0+y2)/h^2
#
#(y1-y2)/(2*h)-evalEQCV(T=5,d=5,sigma=0.5,K=(1-0.0)*expectA(T=5,d=5),S0=100,greeks=TRUE)[2]
#(y1-2*y0+y2)/h^2-evalEQCV(T=5,d=5,sigma=0.5,K=(1-0.0)*expectA(T=5,d=5),S0=100,greeks=TRUE)[3]
#
#h <-1e-3
#y1 <- evalEQCV(T=5,d=5,sigma=0.5,K=(1-0.0)*expectA(T=5,d=5),S0=100+h,greeks=TRUE)[2]
#y2 <- evalEQCV(T=5,d=5,sigma=0.5,K=(1-0.0)*expectA(T=5,d=5),S0=100-h,greeks=TRUE)[2]
#(y1-y2)/(2*h)-evalEQCV(T=5,d=5,sigma=0.5,K=(1-0.0)*expectA(T=5,d=5),S0=100,greeks=TRUE)[3]



AsianCall_NCV_CMC_QCV_greeks <- function(n=10^4,T=1,d=12,K=100,r=0.05,sigma=0.2,S0=100,
				maxiter=100,tol=1e-14,np=100){
	# New simulation algorithm for Asian call options 
	# new CV, conditional Monte Carlo and quadratic CVs
	# Algorithm 7 in Dingec and Hormann (2013) p. 430.
	# n ... number of repetitions for the simulation
	# tol ... tolerance of Newton's method
	# maxiter ...  max number of iterations for Newton's method
	# plotyn ... if TRUE, the scatter plot of y against psi is plotted.
	# np ... sample size of pilot simulation

	n <- n+np
	dt <- T/d
	sumSt <- sumlogSt <-0
	logSt <- log(S0)

	varX <- d*(d+1)*(2*d+1)/6
	mus <- log(S0)+(r-sigma^2/2)*dt*(d+1)/2
	sigmas <- sigma/d*sqrt(dt*varX)
	k <- (log(K)-mus)/sigmas

	# Simulate marginals conditional on Z
	v <- (d:1)/sqrt(varX)
	Z <- matrix(rnorm(n*d),d,n)
	zmatr <- Z - v%o%(v%*%Z)[1,]
	
	a <- sigma*sqrt(dt)*cumsum(v)
	s <- matrix(0,d,n)
	s[1,] <- S0*exp((r-sigma^2/2)*dt+sigma*sqrt(dt)*zmatr[1,])
	for(i in 2:d)
		s[i,] <- s[i-1,]*exp((r-sigma^2/2)*dt+sigma*sqrt(dt)*zmatr[i,])	
	
	# Newton-Raphson method to find the root: b(Z)

	# Function as an argument to Newtons method
	# Returns the function value and its derivative at x	
	fdf <- function(x) {
		mat <- exp(a%o%x)*s
		cbind((1/d)*colSums(mat)-K,(1/d)*colSums(a*mat))
	}
	# Approximate root as a starting point (1st order app)
	x0 <- findrootapp(a,s,k,K,d,n)
	# Newtons method
	nwroot <- newtons_root(x0,fdf,maxiter,tol)
	b <- nwroot[,1]

	# Computing conditional expectation
	m <- matrix(0,d,n)
	for(i in 1:d) m[i,] <- exp(a[i]^2/2)*(pnorm(k-a[i])-pnorm(b-a[i]))
	Int <- (1/d)*colSums(m*s)
	y <- exp(-r*T)*(Int-K*(pnorm(k)-pnorm(b)))
	# delta
	expak <- as.vector(exp(a%o%k))
	Ak <- (1/d)*colSums(expak*s)
	delta <- exp(-r*T)*(Int-(Ak-K)*dnorm(k)/sigmas)/S0
	# gamma
	Apb <- (1/d)*colSums(a*exp(a%o%b)*s)
	Apk <- (1/d)*colSums(a*expak*s)
	C1 <- (K^2*dnorm(b)/Apb-Ak*dnorm(k)/sigmas)/S0^2
	C2 <- (Apk-Ak*sigmas+(Ak-K)*(sigmas-k))*dnorm(k)/(S0*sigmas)^2
	gamma <- exp(-r*T)*(C1 + C2)
	#print(var(C1)/var(C1+C2))
	
	# CV
	bcv <- findbcv(K,T,d,r,sigma,S0)
	mcv <- exp(a^2/2)*(pnorm(k-a)-pnorm(bcv-a))
	Intcv <- (1/d)*colSums(mcv*s)
	cecv <- Intcv-K*(pnorm(k)-pnorm(bcv))
	x1 <- exp(-r*T)* cecv
	x2 <- cecv^2
	# delta CV
	Abcv <- (1/d)*colSums(as.vector(exp(a%o%bcv))*s)
	h0bcv <- (1/d)*sum(exp(a*bcv+r*(1:d)*dt-a^2/2))
	hpbcv <- (S0/d)*sum(a*exp(a*bcv+r*(1:d)*dt-a^2/2))
	bcv0 <- -h0bcv/hpbcv 
	deltacv <- Intcv/S0-(Ak-K)*dnorm(k)/(S0*sigmas)-(Abcv-K)*dnorm(bcv)*bcv0
	#deltacv <- exp(-r*T)*deltacv
	
	# gamma CV
	Apbcv <- (1/d)*colSums(a*as.vector(exp(a%o%bcv))*s)
	hppbcv <- (S0/d)*sum(a^2*exp(a*bcv+r*(1:d)*dt-a^2/2))
	bcv00 <- -bcv0*(2/S0+ hppbcv* bcv0/hpbcv)
	C1 <- -Ak*dnorm(k)/(S0^2*sigmas)-Abcv/S0*dnorm(bcv)*bcv0
	C2 <- dnorm(k)/(S0*sigmas)^2*(Apk-Ak*sigmas+(Ak-K)*(sigmas-k))
	C3 <- dnorm(bcv)*(Abcv/S0*bcv0+ Apbcv*bcv0^2+(Abcv-K)*(bcv00-bcv*bcv0^2))
	gammacv <- (C1+C2-C3)


	# if(plotyn==TRUE){
		# plot(x1,y,cex=.5,xlab=expression(psi),ylab="Y")
		# #lines(sort(cecv),sort(reg$fitted.values))
	# }
	
	# QCVs
	deltacv2 <- 2*deltacv*cecv
	gammacv2 <- 2*(deltacv^2+ cecv*gammacv)
	
	# Expectations
	expectationsLB <- evalLB(K,T,d,r,sigma,S0,greeks=TRUE,all=FALSE)
	ex1 <- expectationsLB[1]
	edeltacv <- exp(r*T)*expectationsLB[2]
	egammacv <- exp(r*T)*expectationsLB[3]
	
	expectationsQCV <- evalEQCV(T,d,K,r,sigma,S0,greeks=TRUE)
	ex2 <- expectationsQCV[1]
	edeltacv2 <- expectationsQCV[2]
	egammacv2 <- expectationsQCV[3]
	
	# Linear Regression
	reg <- lm(cbind(y,delta,gamma)[1:np,]~x1[1:np]+x2[1:np]+deltacv[1:np]+ gammacv[1:np]+ deltacv2[1:np]+ gammacv2[1:np])
	cs <- reg$coeff[-1,]
	
	#print(summary(reg))
	#reg <- lm(y[1:np]~x1[1:np]+x2[1:np])
	#cs <- reg$coeff[-1]
	
	#ex1 <- evalLB(K,T,d,r,sigma,S0,FALSE)
	#ex2 <- eval_equad(T,d,K,r,sigma,S0)
	
	#print(c(mean(deltacv), edeltacv))
	#print(c(mean(deltacv2), edeltacv2))
	#print(c(mean(gammacv), egammacv))
	#print(c(mean(gammacv2), egammacv2))
	
	ycv <- cbind(y,delta,gamma)-cbind(x1-ex1,x2-ex2,deltacv-edeltacv,gammacv-egammacv, deltacv2-edeltacv2, gammacv2-egammacv2)%*%cs

	ew <- evalECV(T,d,K,r,sigma,S0,greeks=TRUE)
	price <- ycv[-(1:np),1] + ew[1]
	delta <- ycv[-(1:np),2] + ew[2]
	gamma <- ycv[-(1:np),3] + ew[3]
	n <- n-np
	
	error_price <- 1.96*sd(price)/sqrt(n)
	error_delta <- 1.96*sd(delta)/sqrt(n)
	error_gamma <- 1.96*sd(gamma)/sqrt(n)
	res <- rbind(c(mean(price),error_price),c(mean(delta),error_delta),c(mean(gamma),error_gamma))
	rownames(res)<-c("price","delta","gamma")
	colnames(res)<-c("result","error estimate")
	res
	
}

#AsianCall_NCV_CMC_QCV_greeks()


#library(randtoolbox)

# Korobov

returnPn_Korobov <- function(n=1021,a=331,d=5){
	# Algorithm of Lemieux 2009, p. 150
	#d <- d+1
	u <- matrix(0,d,n)
	z <- rep(1,d)
	for(j in 2:d) z[j] <- (a*z[j-1])%%n
	for(i in 2:n) u[,i] <-(u[,i-1]+z/n)%%1
	u#[,-1]
}

#returnPn_Korobov(n=1021,a=331,d=5)[,1:2]

# generators

#nvec <- c(1021,4093,16381,65521)
#avec <- c(331,219,665,2469)
#avec <- c(76,1516,4026,8950)

# see Lecuyer and Lemieux (200) for different generators

# naive simulation

simulateAsianCall_Z <- function(Zmat,T=1,K=100,r=0.05,sigma=0.1,S0=100){
	# Naive Simulation Algorithm for Asian option with arithmetic average
	# n ... number of repetitions for the simulation
	d <- ncol(Zmat)
	n <- nrow(Zmat)
	dt <- T/d
	sumSt <- 0
	logSt <- log(S0)
	for(i in 1:d){
		logSt <- logSt+(r-sigma^2/2)*dt+sigma*sqrt(dt)* Zmat[,i]
		sumSt <- sumSt + exp(logSt)
	}
	mean(exp(-r*T)*pmax(sumSt/d-K,0))
}


AsianCall_naive_greeks_Z <- function(Z,T=1,d=12,K=100,r=0.05,sigma=0.1,S0=100,genmethod="pca",Qmat){
	# Naive Simulation Algorithm for Asian option with arithmetic average
	# Algorithm 2 in Dingec and Hormann (2013), p.423
	# n ... number of repetitions for the simulation
	
	d <- nrow(Z)
	n <- ncol(Z)
	
	dt <- T/d
	sumSt <- 0
	
	if(genmethod=="pca"){
		zmatr <- Qmat %*% Z
		St <- S0*exp((r-sigma^2/2)*(1:d)*dt+sigma*sqrt(dt)*zmatr)
		sumSt <- colSums(St)
		Z1 <- zmatr[1,]
	}
	else if (genmethod=="std"){
		logSt <- log(S0)
		for(i in 1:d){
			logSt <- logSt+(r-sigma^2/2)*dt+sigma*sqrt(dt)*Z[i,]
			sumSt <- sumSt + exp(logSt)
		}
		Z1 <- Z[1,]
	}
	else{
		stop("such genmethod does not exist")
	}
	A <- sumSt/d
	price <- exp(-r*T)*pmax(A-K,0)	
	delta <- exp(-r*T)*(A>K)*A/S0
	# mixed estimator (LR-PW) for gamma
	L1 <- sigma*sqrt(dt)
	gamma <- exp(-r*T)*(A>K)*K*Z1/(S0^2*L1)
	y <- cbind(price, delta, gamma)
	colMeans(y)
}

#AsianCall_naive_greeks_Z(matrix(rnorm(12*1000),1000,12))

AsianCall_naive_greeks_qmc <- function(nout=25,n=1024,T=1,d=12,K=100,r=0.05,sigma=0.2,S0=100,scrambling = 2, seed = 4711, genmethod ="pca", dirnum=1, seq.type="sobol", a, baker=TRUE){
	
	y <- matrix(0,3,nout)
	if(seq.type=="korobov") Pn <- returnPn_Korobov(n,a,d)
	
	if(genmethod=="pca") {
		L <- matrix(1,d,d)
		L[upper.tri(L)]<-0
		eigenres <- eigen(L%*%t(L))
		# print(system.time(eigenres <- eigen(L%*%t(L))))
		Qmat <- eigenres$vectors%*%sqrt(diag(eigenres$values))
	}
	
	for(i in 1:nout){
		if(seq.type=="sobol"){
            print("sobol() not activated in the moment")
		#			Z <- t(sobol(n=n, dim = d, scrambling = scrambling, seed = seed, normal = TRUE))
			seed <- seed+1
		}
		else if(seq.type=="korobov"){
			U <- (runif(d)+Pn)%%1
			if(baker){
				Ub <- U
				U[Ub<=0.5] <- 2*U[Ub<=0.5]
				U[Ub>0.5] <- 2*(1-U[Ub>0.5])
			}
			# print(sum(U==1))
			Z <- qnorm(U)
		}
		Z[which(Z==-Inf)] <- qnorm(1/n) # a primitive solution for the problem u=0
		y[,i] <- AsianCall_naive_greeks_Z(Z,T,d=d,K,r,sigma,S0,genmethod,Qmat)
		seed <- seed+1
	}
	error <- 1.96*apply(y,1,sd)/sqrt(nout)
	res <- cbind(rowMeans(y), error)
	colnames(res)<- c("result","error estimate")
	rownames(res) <- c("price","delta","gamma")
	res
}

#AsianCall_naive_greeks_qmc()
#AsianCall_naive_greeks(n=1024*25)
#(AsianCall_naive_greeks(n=1024*25,sigma=0.2)/AsianCall_naive_greeks_qmc())^2


AsianCall_NCV_CMC_greeks_Z <- function(Z,T=1,d=12,K=100,r=0.05,sigma=0.1,S0=100,
				maxiter=100,tol=1e-14,pca=FALSE){
	# New simulation algorithm for Asian call options 
	# new CV, conditional Monte Carlo and quadratic CVs
	# Algorithm 7 in Dingec and Hormann (2013) p. 430.
	# n ... number of repetitions for the simulation
	# tol ... tolerance of Newton's method
	# maxiter ...  max number of iterations for Newton's method
	# plotyn ... if TRUE, the scatter plot of y against psi is plotted.
	# np ... sample size of pilot simulation

	d <- nrow(Z)
	n <- ncol(Z)

	dt <- T/d
	sumSt <- sumlogSt <-0
	logSt <- log(S0)

	varX <- d*(d+1)*(2*d+1)/6
	mus <- log(S0)+(r-sigma^2/2)*dt*(d+1)/2
	sigmas <- sigma/d*sqrt(dt*varX)
	k <- (log(K)-mus)/sigmas

	# Simulate marginals conditional on Z
	v <- (d:1)/sqrt(varX)	
	a <- sigma*sqrt(dt)*cumsum(v)

	if(pca){
		L <- matrix(1,d,d)
		L[upper.tri(L)]<-0
		A <- L%*%(diag(1,d)-v%o%v)
		eigenres <- eigen(A%*%t(A),symmetric=TRUE)
		V <- eigenres$vectors
		lam <- eigenres$values
		lam[d]<-0
		Lmbd <- diag(lam)
		A <- V%*%sqrt(Lmbd)
		zmatr <- A%*%Z	
		s <- S0*exp((r-sigma^2/2)*(1:d)*dt+sigma*sqrt(dt)*zmatr)
	}
	else{
		zmatr <- Z - v%o%(v%*%Z)[1,]
		s <- matrix(0,d,n)
		s[1,] <- S0*exp((r-sigma^2/2)*dt+sigma*sqrt(dt)*zmatr[1,])
		for(i in 2:d)
			s[i,] <- s[i-1,]*exp((r-sigma^2/2)*dt+sigma*sqrt(dt)*zmatr[i,])	
	}
	
	# Newton-Raphson method to find the root: b(Z)

	# Function as an argument to Newtons method
	# Returns the function value and its derivative at x	
	fdf <- function(x) {
		mat <- exp(a%o%x)*s
		cbind((1/d)*colSums(mat)-K,(1/d)*colSums(a*mat))
	}
	# Approximate root as a starting point (1st order app)
	x0 <- findrootapp(a,s,k,K,d,n)
	if(is.nan(max(x0))) print(which(is.nan(x0)))
	# Newtons method
	nwroot <- newtons_root(x0,fdf,maxiter,tol)
	b <- nwroot[,1]

	# Computing conditional expectation
	m <- matrix(0,d,n)
	for(i in 1:d) m[i,] <- exp(a[i]^2/2)*(pnorm(k-a[i])-pnorm(b-a[i]))
	Int <- (1/d)*colSums(m*s)
	y <- exp(-r*T)*(Int-K*(pnorm(k)-pnorm(b)))
	# delta
	expak <- as.vector(exp(a%o%k))
	Ak <- (1/d)*colSums(expak*s)
	delta <- exp(-r*T)*(Int-(Ak-K)*dnorm(k)/sigmas)/S0
	# gamma
	Apb <- (1/d)*colSums(a*exp(a%o%b)*s)
	Apk <- (1/d)*colSums(a*expak*s)
	C1 <- (K^2*dnorm(b)/Apb-Ak*dnorm(k)/sigmas)/S0^2
	#C1 <- (K^2*dnorm(b)/Apb)/S0^2
	C2 <- (Apk-Ak*sigmas+(Ak-K)*(sigmas-k))*dnorm(k)/(S0*sigmas)^2
	gamma <- exp(-r*T)*(C1 + C2)
	
	# CV
	bcv <- findbcv(K,T,d,r,sigma,S0)
	mcv <- exp(a^2/2)*(pnorm(k-a)-pnorm(bcv-a))
	Intcv <- (1/d)*colSums(mcv*s)
	cecv <- Intcv-K*(pnorm(k)-pnorm(bcv))
	x1 <- exp(-r*T)* cecv
	x2 <- cecv^2

	y <- rbind(y, delta, gamma)
	res1 <- rowMeans(y)
	
	# delta CV
	Abcv <- (1/d)*colSums(as.vector(exp(a%o%bcv))*s)
	h0bcv <- (1/d)*sum(exp(a*bcv+r*(1:d)*dt-a^2/2))
	hpbcv <- (S0/d)*sum(a*exp(a*bcv+r*(1:d)*dt-a^2/2))
	bcv0 <- -h0bcv/hpbcv 
	deltacv <- Intcv/S0-(Ak-K)*dnorm(k)/(S0*sigmas)-(Abcv-K)*dnorm(bcv)*bcv0
	deltacv <- exp(-r*T)*deltacv

	# gamma CV
	Apbcv <- (1/d)*colSums(a*as.vector(exp(a%o%bcv))*s)
	h0pbcv <- hpbcv/S0
	bcv00 <- h0bcv*h0pbcv/hpbcv^2
	C1 <- -Ak*dnorm(k)/(S0^2*sigmas)-Abcv/S0*dnorm(bcv)*bcv0
	C2 <- dnorm(k)/(S0*sigmas)^2*(Apk-Ak*sigmas+(Ak-K)*(sigmas-k))
	C3 <- dnorm(bcv)*(Abcv/S0*bcv0+ Apbcv*bcv0^2+(Ak-K)*(bcv00-bcv*bcv0^2))
	gammacv <- exp(-r*T)*(C1+C2-C3)

	x <- rbind(x1, x2, deltacv, gammacv, deltacv*x1, deltacv^2+x1*gammacv)
	res2<- rowMeans(x)
	list(res1,res2)
}

AsianCall_NCV_CMC_greeks_qmc <- function(nout=25,noutp=25,n=1024,T=1,d=12,K=100,r=0.05,sigma=0.1,S0=100,maxiter=100,tol=1e-14,scrambling = 2, seed = 4711,pca=TRUE,seq.type="sobol",a){
	
	# Pilot run
	y <- matrix(0,3,noutp)
	x1p <- x2p <- rep(0,noutp)
	deltacv2 <- gammacv2 <- deltacv <- gammacv <- rep(0,noutp)
	# For korobov
	if(seq.type=="korobov") Pn <- returnPn_Korobov(n,a,d)
	for(i in 1:noutp){
		if(seq.type=="sobol"){
            print("sobol() not activated in the moment")
		#				Z <- t(sobol(n=n, dim = d, scrambling = scrambling, seed = seed, normal = TRUE))
			seed <- seed+1
		}
		else if(seq.type=="korobov"){
			Z <- qnorm((runif(d)+Pn)%%1)
		}
		Z[which(Z==-Inf)] <- qnorm(1/n) # a primitive solution for the problem u=0
		res <- AsianCall_NCV_CMC_greeks_Z(Z,T,d,K,r,sigma,S0,maxiter,tol,pca)
		y[,i] <- res[[1]]
		x1p[i] <- res[[2]][1]
		x2p[i] <- res[[2]][2]
		deltacv[i] <- res[[2]][3]  
		gammacv[i] <- res[[2]][4]
		deltacv2[i] <- res[[2]][5]
		gammacv2[i] <- res[[2]][6]  
		#print(c(i,seed))
	}
	reg <- lm(y[1,]~x1p+x2p)
	csp <- reg$coefficients[-1]
	#print("price")
	#print(1/(1-summary(reg)$r.squared))
	#print(1/(1-summary(lm(y[1,]~x1p+x2p+ deltacv+ gammacv+ deltacv2+ gammacv2))$r.squared))
	print((1/(1-summary(lm(y[1,]~x1p+x2p+ deltacv+ gammacv+ deltacv2+ gammacv2))$r.squared))/(1/(1-summary(reg)$r.squared)))

	#plot(deltacv2,y[2,])
	#print(1/(1-cor(y[2,], deltacv2)^2))
	#print("delta")
	#print(1/(1-summary(lm(y[2,]~ deltacv+ deltacv2))$r.squared))
	#print(1/(1-summary(lm(y[2,]~ deltacv+ deltacv2+gammacv+ gammacv2))$r.squared))
	print(1/(1-summary(lm(y[2,]~ deltacv+ deltacv2+gammacv+ gammacv2+x1p+x2p))$r.squared))
	#print("gamma")
	#print(1/(1-summary(lm(y[3,]~ gammacv+ gammacv2))$r.squared))
	#print(1/(1-summary(lm(y[3,]~ gammacv+ gammacv2+deltacv+ deltacv2))$r.squared))
	print(1/(1-summary(lm(y[3,]~ gammacv+ gammacv2+deltacv+ deltacv2 +x1p+x2p))$r.squared))


	# Expectations
	ew <- evalECV(T,d,K,r,sigma,S0,TRUE)
	ex1 <- evalLB(K,T,d,r,sigma,S0,FALSE)
	ex2 <- eval_equad(T,d,K,r,sigma,S0)

	
	# Main simulation
	y <- matrix(0,3,nout)
	x1p <- x2p <-rep(0,nout)
	for(i in 1:nout){
		if(seq.type=="sobol"){
            print("sobol() not activated in the moment")
		#				Z <- t(sobol(n=n, dim = d, scrambling = scrambling, seed = seed, normal = TRUE))
			seed <- seed+1
		}
		else if(seq.type=="korobov"){
			Z <- qnorm((runif(d)+Pn)%%1)
		}
		Z[which(Z==-Inf)] <- qnorm(1/n) # a primtive solution for the problem u=0
		res <- AsianCall_NCV_CMC_greeks_Z(Z,T,d,K,r,sigma,S0,maxiter,tol,pca)
		y[,i] <- res[[1]]
		x1p[i] <- res[[2]][1]
		x2p[i] <- res[[2]][2]
	}
	y <- y +ew
	#y[-3,] <- y[-3,] +ew[-3]
	y[1,] <- y[1,] -csp[1]*(x1p-ex1)- csp[2]*(x2p-ex2)
	error <- 1.96*apply(y,1,sd)/sqrt(nout)
	res <- cbind(rowMeans(y), error)
	colnames(res)<- c("result","error estimate")
	rownames(res) <- c("price","delta","gamma")
	res
}

#AsianCall_NCV_CMC_greeks_qmc(seq.type="sobol",n=1024)
#AsianCall_NCV_CMC_greeks_qmc(seq.type="korobov",n=1021,a=76)


#AsianCall_CMC_CV_greeks_qmc(n=2^16,scrambling=2,seed=4721)





# problem: getting u=0 for sramb=2 and large n

# sobol(n=2^16,dim=12,scrambling=2,seed=4720)[36739,]
# sobol(n=2^16,dim=12,scrambling=2,seed=4745)[61591,]



# Combination with QCV


AsianCall_NCV_CMC_QCV_greeks_Z <- function(Z,T=1,d=12,K=100,r=0.05,sigma=0.1,S0=100,
				maxiter=100,tol=1e-14,genmethod="pca",Qmat){
	# New simulation algorithm for Asian call options 
	# new CV, conditional Monte Carlo and quadratic CVs
	# Algorithm 7 in Dingec and Hormann (2013) p. 430.
	# n ... number of repetitions for the simulation
	# tol ... tolerance of Newton's method
	# maxiter ...  max number of iterations for Newton's method
	# plotyn ... if TRUE, the scatter plot of y against psi is plotted.
	# np ... sample size of pilot simulation

	d <- nrow(Z)
	n <- ncol(Z)

	dt <- T/d
	varX <- d*(d+1)*(2*d+1)/6
	mus <- log(S0)+(r-sigma^2/2)*dt*(d+1)/2
	sigmas <- sigma/d*sqrt(dt*varX)
	k <- (log(K)-mus)/sigmas

	# Simulate marginals conditional on Z
	v <- (d:1)/sqrt(varX)	
	a <- sigma*sqrt(dt)*cumsum(v)
	
	if(genmethod=="std"){
		zmatr <- Z - v%o%(v%*%Z)[1,]
		s <- matrix(0,d,n)
		s[1,] <- S0*exp((r-sigma^2/2)*dt+sigma*sqrt(dt)*zmatr[1,])
		for(i in 2:d)
			s[i,] <- s[i-1,]*exp((r-sigma^2/2)*dt+sigma*sqrt(dt)*zmatr[i,])	
	}
	else {
		if(genmethod=="pca"||genmethod=="ltpca") Z<-Z[-d,]
		#print(dim(Qmat))
		#print(dim(Z))
		zmatr <- Qmat%*%Z	
		s <- S0*exp((r-sigma^2/2)*(1:d)*dt+sigma*sqrt(dt)*zmatr)	
	}

	
	# Newton-Raphson method to find the root: b(Z)

	# Function as an argument to Newtons method
	# Returns the function value and its derivative at x	
	fdf <- function(x) {
		mat <- exp(a%o%x)*s
		cbind((1/d)*colSums(mat)-K,(1/d)*colSums(a*mat))
	}
	# Approximate root as a starting point (1st order app)
	x0 <- findrootapp(a,s,k,K,d,n)
	if(is.nan(max(x0))) {
		print(which(is.nan(x0)))
	}
	# Newtons method
	nwroot <- newtons_root(x0,fdf,maxiter,tol)
	b <- nwroot[,1]

	# Computing conditional expectation
	m <- matrix(0,d,n)
	for(i in 1:d) m[i,] <- exp(a[i]^2/2)*(pnorm(k-a[i])-pnorm(b-a[i]))
	Int <- (1/d)*colSums(m*s)
	y <- exp(-r*T)*(Int-K*(pnorm(k)-pnorm(b)))
	# delta
	expak <- as.vector(exp(a%o%k))
	Ak <- (1/d)*colSums(expak*s)
	delta <- exp(-r*T)*(Int-(Ak-K)*dnorm(k)/sigmas)/S0
	# gamma
	Apb <- (1/d)*colSums(a*exp(a%o%b)*s)
	Apk <- (1/d)*colSums(a*expak*s)
	C1 <- (K^2*dnorm(b)/Apb-Ak*dnorm(k)/sigmas)/S0^2
	#C1 <- (K^2*dnorm(b)/Apb)/S0^2
	C2 <- (Apk-Ak*sigmas+(Ak-K)*(sigmas-k))*dnorm(k)/(S0*sigmas)^2
	gamma <- exp(-r*T)*(C1 + C2)
	
	# CV
	bcv <- findbcv(K,T,d,r,sigma,S0)
	mcv <- exp(a^2/2)*(pnorm(k-a)-pnorm(bcv-a))
	Intcv <- (1/d)*colSums(mcv*s)
	cecv <- Intcv-K*(pnorm(k)-pnorm(bcv))
	x1 <- exp(-r*T)* cecv
	x2 <- cecv^2
	
	# delta CV
	Abcv <- (1/d)*colSums(as.vector(exp(a%o%bcv))*s)
	h0bcv <- (1/d)*sum(exp(a*bcv+r*(1:d)*dt-a^2/2))
	hpbcv <- (S0/d)*sum(a*exp(a*bcv+r*(1:d)*dt-a^2/2))
	bcv0 <- -h0bcv/hpbcv 
	deltacv <- Intcv/S0-(Ak-K)*dnorm(k)/(S0*sigmas)-(Abcv-K)*dnorm(bcv)*bcv0
	#deltacv <- exp(-r*T)*deltacv
	
	# gamma CV
	Apbcv <- (1/d)*colSums(a*as.vector(exp(a%o%bcv))*s)
	hppbcv <- (S0/d)*sum(a^2*exp(a*bcv+r*(1:d)*dt-a^2/2))
	bcv00 <- -bcv0*(2/S0+ hppbcv* bcv0/hpbcv)
	C1 <- -Ak*dnorm(k)/(S0^2*sigmas)-Abcv/S0*dnorm(bcv)*bcv0
	C2 <- dnorm(k)/(S0*sigmas)^2*(Apk-Ak*sigmas+(Ak-K)*(sigmas-k))
	C3 <- dnorm(bcv)*(Abcv/S0*bcv0+ Apbcv*bcv0^2+(Abcv-K)*(bcv00-bcv*bcv0^2))
	gammacv <- C1+C2-C3
	
	deltacv2 <- 2*deltacv*cecv
	gammacv2 <- 2*(deltacv^2+ cecv*gammacv)

	y <- rbind(y, delta, gamma)
	res1 <- rowMeans(y)
	x <- rbind(x1, x2, deltacv, gammacv, deltacv2, gammacv2)
	res2 <- rowMeans(x)
	list(res1,res2)
}




ortmat<-function(mat)
{
  n<-dim(mat)[1]
  k<-dim(mat)[2]
  res<-matrix(0,n,n)
  res[1:n,1:k]<-mat
  
  for(i in (k+1):n){
    res[i:n,i]<-1
    A<-t(res[1:(i-1),1:(i-1)])
    b<-t(-(t(res[i:n,i])%*%res[i:n,1:(i-1)]))
    res[1:(i-1),i]<-solve(A,b)
  }
  for(i in 1:n){
    res[,i]<-res[,i]/sqrt(sum(res[,i]^2))
  }
  res
}


findQmat <- function(T=1,d=12,K=100,r=0.05,sigma=0.2,S0=100, genmethod="pca", dirnum=1){
	dt <- T/d
	varX <- d*(d+1)*(2*d+1)/6
	mus <- log(S0)+(r-sigma^2/2)*dt*(d+1)/2
	sigmas <- sigma/d*sqrt(dt*varX)
	k <- (log(K)-mus)/sigmas
	v <- (d:1)/sqrt(varX)	
	a <- sigma*sqrt(dt)*cumsum(v)
	L <- matrix(1,d,d)
	L[upper.tri(L)]<-0
	B <- L%*%(diag(1,d)-v%o%v)
	# generation methods
	if(genmethod=="pca"){
		eigenres <- eigen(B%*%t(B),symmetric=TRUE)		
		V <- eigenres$vectors
		lam <- eigenres$values
		lam[d]<-0
		Lmbd <- diag(lam)
		Qmat <- (V%*%sqrt(Lmbd))[,-d]
	}
	else if(genmethod=="pcamain"){
		# using only first few main directions
		eigenres <- eigen(B%*%t(B),symmetric=TRUE)	
		V <- eigenres$vectors
		v <- t(B) %*% V[, 1:dirnum]
		# A <- matrix(1,d,d)
		# A[,1:dirnum] <- v
		# A <- qr.Q(qr(A))
		A <- ortmat(v)
		Qmat <- B%*%A
	}
	else if (genmethod=="lt"){
		eps <- rep(0,d)
		A <- matrix(1,d,d)
		for(i in 1:dirnum){
			if(i>1) eps[1:(i-1)] <- 1
			s <- S0*exp((r-sigma^2/2)*(1:d)*dt+sigma*sqrt(dt)*as.vector(B%*%A%*%eps))
			b <- uniroot(function(x) 1/d*sum(exp(a*x)*s)-K,c(-10,k))$root
			m <- exp(a^2/2)*(pnorm(k-a)-pnorm(b-a))
			dh <- sigma*sqrt(dt)/d*colSums(m*s*B)
			A[,i] <- dh
			A[,1:i] <- qr.Q(qr(A[,1:i]))
		}
		if(dirnum<d)  {
			#A <- qr.Q(qr(A))
			A <- ortmat(A[,1: dirnum])
		}
		Qmat <- B%*%A
	}
	else if (genmethod=="ltpca"){
		# combination of lt with pca
		eigenres <- eigen(B%*%t(B),symmetric=TRUE)
		V <- eigenres$vectors
		lam <- eigenres$values
		lam[d]<-0
		Lmbd <- diag(lam)
		B <- (V%*%sqrt(Lmbd))[,-d]
		eps <- rep(0,d-1)
		A <- matrix(1,d-1,d-1)
		for(i in 1:(d-1)){
			if(i>1) eps[1:(i-1)] <- 1
			s <- S0*exp((r-sigma^2/2)*(1:d)*dt+sigma*sqrt(dt)*as.vector(B%*%A%*%eps))
			b <- uniroot(function(x) 1/d*sum(exp(a*x)*s)-K,c(-10,k))$root
			m <- exp(a^2/2)*(pnorm(k-a)-pnorm(b-a))
			dh <- sigma*sqrt(dt)/d*colSums(m*s*B)
			A[,i] <- dh
			A[,1:i] <- qr.Q(qr(A[,1:i]))
		}
		Qmat <- B%*%A	
	}
	Qmat
}

#system.time(findQmat(T=1,d=12,K=100,r=0.05,sigma=0.2,S0=100, genmethod="pca", dirnum=1))


AsianCall_NCV_CMC_QCV_greeks_qmc <- function(nout=25,noutp=25,n=1024,T=1,d=12,K=100,r=0.05,sigma=0.2,S0=100,maxiter=100,tol=1e-14, scrambling =2, seed = 4711, genmethod="pca", dirnum=1, seq.type="sobol", a, baker=TRUE, cvmethod="splitting", makehist=FALSE){
	
	pilotyn <- cvmethod=="pilotrun"
	y <- matrix(0, pilotyn*noutp+nout,3)
	X <- matrix(0, pilotyn*noutp+nout,6)
	if(seq.type=="korobov") Pn <- returnPn_Korobov(n,a,d)
	
	if(genmethod!="std") Qmat <- findQmat(T,d,K,r,sigma,S0,genmethod,dirnum)

	# Outer repetitions
	
	for(i in 1:(pilotyn*noutp+nout)){
		if(seq.type=="sobol"){
            print("sobol() not activated in the moment")
		#				Z <- t(sobol(n=n, dim = d, scrambling = scrambling, seed = seed, normal = TRUE))
			seed <- seed+1
		}
		else if(seq.type=="korobov"){
			U <- (runif(d)+Pn)%%1
			if(baker){
				Ub <- U
				U[Ub<=0.5] <- 2*U[Ub<=0.5]
				U[Ub>0.5] <- 2*(1-U[Ub>0.5])
			}
			# print(sum(U==1))
			Z <- qnorm(U)
		}
		Z[which(Z==-Inf)] <- qnorm(1/n) # a primitive solution for the problem u=0
		Z[which(Z==Inf)] <- qnorm(1-1/n) # a primitive solution for the problem u=1 
		res <- AsianCall_NCV_CMC_QCV_greeks_Z(Z,T,d,K,r,sigma,S0,maxiter,tol,genmethod,Qmat)
		# print(res)
		y[i,] <- res[[1]]
		X[i,] <- res[[2]]
		#print(c(i,seed))
	}
	# plot(y[-1,1],y[-nout,1]) # shows independence

	# Expectations
	ew <- evalECV(T,d,K,r,sigma,S0,greeks=TRUE)
	expectationsLB <- evalLB(K,T,d,r,sigma,S0,greeks=TRUE,all=FALSE)
	ex1 <- expectationsLB[1]
	edeltacv <- exp(r*T)*expectationsLB[2]
	egammacv <- exp(r*T)*expectationsLB[3]
	
	expectationsQCV <- evalEQCV(T,d,K,r,sigma,S0,greeks=TRUE)
	ex2 <- expectationsQCV[1]
	edeltacv2 <- expectationsQCV[2]
	egammacv2 <- expectationsQCV[3]
	
	EX <- c(ex1, ex2, edeltacv, egammacv, edeltacv2, egammacv2)
	
	# Linear Regression
	if(cvmethod=="direct"){
		cs <- lm(y~X)$coeff[-1,]
		# print(summary(lm(y~X)))
		# print(cs)
		# print(cbind(lm(y[,1]~X)$coeff[-1],lm(y[,2]~X)$coeff[-1],lm(y[,3]~X)$coeff[-1])==cs)
		ycv <- y - t(t(X)-EX)%*% cs
	}
	if(cvmethod=="splitting"){
		ycv <- matrix(0,nout,3)
		#print(system.time(
		for (i in 1:nout){
			reg <- lm(y[-i,]~X[-i,])
			cs <- reg$coeff[-1,]
			ycv[i,] <- y[i,] - (X[i,]-EX)%*% cs
		}
		#))
		# bias estimation of the direct method
		# cs <- lm(y~X)$coeff[-1,]
		# ycv <- ycv - (y - t(t(X)-EX)%*% cs)
	}
	else if (cvmethod=="pilotrun"){
		reg <- lm(y[1:noutp,]~X[1:noutp,])
		cs <- reg$coeff[-1,]
		ycv <- y - t(t(X)-EX)%*% cs
		ycv <- ycv[-(1:noutp),]
	}
	#print(cs)
	ycv <-  t(ycv) + ew
	# plot(ycv[2,-1],ycv[2,-nout]) # shows independence
	if(makehist){
		hist(ycv[1,],breaks=50)
	}
	error <- 1.96*apply(ycv,1,sd)/sqrt(nout)
	res <- cbind(rowMeans(ycv), error)
	colnames(res)<- c("result","error estimate")
	rownames(res) <- c("price","delta","gamma")
	res
}

#AsianCall_NCV_CMC_QCV_greeks_qmc()

#############################
#
# wrapper function that is exported
#AsianCall_naive_greeks <- function(n=10^4,T=1,d=12,K=100,r=0.05,sigma=0.2,S0=100)
#AsianCall_naive_greeks_qmc <- function(nout=25,n=1024,T=1,d=12,K=100,r=0.05,sigma=0.2,S0=100,scrambling = 2, seed = 4711, genmethod ="pca", 
#                             dirnum=1, seq.type="sobol", a, baker=TRUE)
#AsianCall_NCV_CMC_QCV_greeks <- function(n=10^4,T=1,d=12,K=100,r=0.05,sigma=0.2,S0=100,maxiter=100,tol=1e-14,np=100)
#AsianCall_NCV_CMC_QCV_greeks_qmc <- function(nout=25,noutp=25,n=1024,T=1,d=12,K=100,r=0.05,sigma=0.2,S0=100,maxiter=100,tol=1e-14, scrambling =2, seed = 4711, 
#                                     genmethod="pca", dirnum=1, seq.type="sobol", a, baker=TRUE, cvmethod="splitting", makehist=FALSE)

AsianCall <- function(T=1,d=12,K=100,r=0.05,sigma=0.2,S0=100,method=c("best","naive"),sampling=c("QMC","MC"),
                       metpar=list(maxiter=100,tol=1.e-14,cvmethod="splitting"),
                       sampar=list(nout=50,seq.type="korobov",n=2039,a=1487,baker=TRUE,genmethod="pca")){
if(method[1]=="best"&sampling[1]=="QMC"){
 res <- AsianCall_NCV_CMC_QCV_greeks_qmc(nout=sampar$nout,n=sampar$n,T=T,d=d,K=K,r=r,sigma=sigma,S0=S0,maxiter=metpar$maxiter,
                                         tol=metpar$tol,genmethod=sampar$genmethod,dirnum=sampar$dirnum,seq.type="korobov",
										 a=sampar$a, baker=sampar$baker,cvmethod=metpar$cvmethod)
}else if(method[1]=="best"&sampling[1]=="MC"){
 res <- AsianCall_NCV_CMC_QCV_greeks(n=sampar$n,T=T,d=d,K=K,r=r,sigma=sigma,S0=S0,maxiter=metpar$maxiter,tol=metpar$tol,
                                     np=metpar$np)
}else if(method[1]=="naive"&sampling[1]=="QMC"){
 res <- AsianCall_naive_greeks_qmc(nout=sampar$nout,n=sampar$n,T=T,d=d,K=K,r=r,sigma=sigma,S0=S0,
                                   genmethod=sampar$genmethod,dirnum=sampar$dirnum,seq.type="korobov",
								   a=sampar$a, baker=sampar$baker)
}else if(method[1]=="naive"&sampling[1]=="MC"){
 res <- AsianCall_naive_greeks(n=sampar$n,T=T,d=d,K=K,r=r,sigma=sigma,S0=S0)
}
return(res)					   
					   
# T=1,d=12,K=100,r=0.05,sigma=0.2,S0=100
# method=c("best","naive")									 
# sampling=c("QMC","MC")
# metpar=list()
#  for method="best"
#   maxiter ... maximal no of iterations for Newton method
#   tol=1.e-14 ... error tolerance for Newton method
#   for sampling="QMC"  
#    cvmethod  ... c("splitting","direct")  NOT necessary for method = "naive"
#  		  "splitting" ... estimates CV coefficients using lm with bootstrap
#         "direct"  ... estimates CV coefficients using lm and the full sample
#   for sampling="MC"  
#      np ... sample size for pilot run for CV; NOT necessary for method = "naive"
# sampar=list()  
#   for sampling="MC" 
#     n ... total samplesize
#   for sampling="QMC" sampar elements:
#     nout  ... number of independent "randomized" copies of the QMC sequence (for QMC)
#     n  ... length of the QMC sequence (forQMC)
#     a ... important constant for the Korobov lattice construction
#     baker ... TRUE/FALSE, indicates if Baker transform should be used for making the integrand periodic
#     genmethod ... c("pca", "std","pcamain","lt",ltpca")
#       note that for method=="naive" only genmethod=c("pca","std") can be used
#       "pca" ... principal component analysis
#       "std" ... standard
#       "pcamain" ... use only first "dirnum" main directions of the PCA
#       "lt" ... uses a transform for the first dirnum
#       "ltpca" ... combination of lt with pca
#        dirnum ... number of main directions, only used for genmethod="pcamain" or "lt"
}

# constants for korobov
#nvec <- c(1021,2039,4093,8191,16381,32749,65521)
#avec1 <- c(331,393,219,1716,665,9515,2469) # worst
#avec2 <- c(76,1487,1516,5130,4026,14251,8950) # optimal
#avec3 <- c(306,280,1397,7151,5693,8363,944)
# avec is from Lecuyer and Lemieux (200)












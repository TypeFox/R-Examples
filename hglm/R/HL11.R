#Updated by Lars 2014 May 28 and 30
HL11 <- function(fv, w, Z, family, tau) {
	n <- nrow(Z)
	p <- ncol(Z)
	mu <- fv
	eta <- family$linkfun(mu)
	dmu.deta <- family$mu.eta(eta)
	g.prime <- 1/dmu.deta
	V <- family$variance(mu)
	#W1 <- (1/g.prime)^2*(1/V) #Bug corrected May 2014
	W1 <- diag((1/g.prime)^2*(1/V)/tau)
	W2 <- diag(w[-(1:n)]^2)
	invC <- solve(t(Z)%*%W1%*%Z+W2)
	K <- Z%*%invC%*%t(Z)
	P <- diag(K)%*%W1
	k1 <- -K
	W1 <- diag(W1)
	#Need to specify the first derivative of the variance function for each possible distribution
	if (family$family[1] == "gaussian" ) V.prime <- rep(0, length(mu))
	if (family$family[1] == "poisson" ) V.prime <- rep(1, length(mu))
	if (family$family[1] == "Gamma" ) V.prime <- 2*mu #Bug corrected May 2014
	if (family$family[1] == "binomial" ) V.prime <- 1 - 2*mu
	#Need to specify the second derivative of the link function for each possible link
	if (family$link == "log")  g.bis <- -1/(mu^2)
	if (family$link == "identity") g.bis <- rep(0,length(mu))
	if (family$link == "inverse") g.bis <- 2/(mu^3)
	if (family$link == "logit") g.bis <- -(1 - 2*mu)*(g.prime^2)
	if (family$link == "probit") g.bis <- eta*(g.prime^2)
	if (family$link == "cloglog") g.bis <- -(log(1 - mu) + 1)*(g.prime^2)
	#dW1.dmu <- -2/(g.prime^3)*g.bis/V - 1/(V^2)*V.prime*(1/g.prime^2) #Bug corrected May 2014
	dW1.dmu <- ( -2/(g.prime^3)*g.bis/V - 1/(V^2)*V.prime*(1/g.prime^2) )/tau
	di <- P[1:n]*(1/W1)^2*dW1.dmu*dmu.deta + crossprod(P[1:n], ((1/W1)*dW1.dmu*dmu.deta*k1[1:n,]))
	return(as.numeric(di/2))
}


#Second derivatives of the link function, ie g.bis, derived by using the following rules.
#Define dmu.deta=f(eta), ie the first derivative of mu as a function of eta
#Then g.bis is: 
#  -[1/(f(eta)^3)]*f'(eta)
#where f'(eta) is the second derivative of mu as a function of eta
#Furthermore, g.prime=1/f(eta)
#####
#For the probit model we have:
# f(eta)=1/sqrt(2*pi)*exp(-0.5*eta^2)
# f'(eta)= -eta*f(eta)
#####
# The derivatives of the logit model can be derived directly since g(mu)=log(mu/(1-mu)) 
# And for the cloglog model we have g(mu)=log(-log(1-mu))

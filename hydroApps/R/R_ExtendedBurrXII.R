
####################################################################################################
# Standard equations of the Extended Burr type XII

# Probability distribution
pBurrXII <- function(x, lambda, k, c) 1-(1-k*(x/lambda)^c)^(1/k)
# Density distribution
dBurrXII <- function(x, lambda, k, c) (c*(1-k*(x/lambda)^c)^(-1+1/k)*(x/lambda)^(-1+c))/lambda
# Quantile function
qBurrXII <- function(p, lambda, k, c) lambda*( 1/k*(1-(1-p)^k) )^(1/c)


####################################################################################################
# Other relationships mostly used as internal functions, but useful for static plots and data interpretation

# Case k->0 (Weibull distribution), lower bound of the domain of non-positive k (tau3 as a function of tau) 
tau3BurrXII.WeibullBound <- function(tau) (1/tau*(-2+2*3^(log(1-tau)/log(2))+3*tau))
# Case k->-inf (Pareto distribution), upper bound of the domain of non-positive k (tau3 as a function of tau) 
tau3BurrXII.ParetoBound <- function(tau) ((1+3*tau)/(3+tau))


####################################################################################################
# Lmoments (i.e. L1=mean, tau=L-CV=L2/L1, tau3=L-skewness=L3/L2) as a function of the distribution parameters
lmomBurrXII <- function(lambda, k, c){
	
	if ((c < 0) | (c < -k)) stop("c must be >0 and >-k")
	if (lambda <= 0) stop("lambda must be >0")
	
	# L1 (mean)
	if (k < 0) {
		L1 <- lambda * (-k)^(-1-(1/c)) * ( gamma(1+1/c)*gamma(-1/c-1/k) ) / ( gamma(1-1/k) )
	} else if (k > 0) {
		L1 <- lambda * k^(-1-(1/c)) * ( gamma(1+1/c)*gamma(1/k) ) / ( gamma(1+1/c+1/k) )
	} else {
		L1 <- lambda*gamma(1+1/c)
	}
	
	# tau = L2/L1
	if (k < 0) {
		tau <- 1 - ( 2*gamma(-(1/c)-2/k) * gamma((-1+k)/k) ) / ( gamma((-2+k)/k) * gamma(-((c+k)/(c*k))) )	
	} else if (k > 0) {
		tau <- 1 - ( 2*gamma(1+1/c+1/k) * gamma(2/k) ) / (gamma(1+1/c+2/k) * gamma(1/k))	
	} else {
		tau <- 1 - 2^(-1/c)
	}
	
	# tau_3 = L3/L2
	if (k < 0) {
		tau3 <- ( (2*gamma(-(1/c)-3/k))/gamma(-(3/k)) - (3*gamma(-(1/c)-2/k))/gamma(-(2/k)) + gamma(-(1/c)-1/k)/gamma(-(1/k)) ) / ( -(gamma(-(1/c)-2/k)/gamma(-(2/k))) + gamma(-(1/c)-1/k)/gamma(-(1/k)) )
	} else if (k > 0) {
		tau3 <- ( gamma(1/k)/gamma(1+1/c+1/k) - (6*gamma(2/k))/gamma(1+1/c+2/k) + (6*gamma(3/k))/gamma(1+1/c+3/k) ) / ( gamma(1/k)/gamma(1+1/c+1/k) - (2*gamma(2/k))/gamma(1+1/c+2/k) )
	} else {
		tau3 <- -( -3 + 2^(1/c) + 2^(1+1/c)*3^(-1/c) ) / ( 1 - 2^(1/c) )
	}

return(c(L1=L1, tau=tau, tau3=tau3))	
	
}

####################################################################################################
# Parameters as a function of the Lmoments. Parameters k and c are appoximated accordi to Ganora & Laio (2014).
# These equations are valid only for pairs of tau-tau3 in the k<=0 domain

parBurrXII.approx <- function(L1, tau, tau3){
	
	if ((tau3 > tau3BurrXII.ParetoBound(tau)) | (tau3 < tau3BurrXII.WeibullBound(tau)))	stop("tau-tau3 outside the validity domain of k<=0")

	# estimate of k (from approximated equation)
	k <- log( -(24-33*tau+9*tau^3-8*(1-tau)^(log(3)/log(2))*(3+tau)+sqrt((24-33*tau+9*tau^3-8*(1-tau)^(log(3)/log(2))*(3+tau))^2-24*tau*(-4*(1-tau)^(log(3)/log(2))*(3+tau)+3*(4-5*tau+tau^3))*(-1+tau*(-3+tau3)+3*tau3)))/(24*(1-tau)^(log(3)/log(2))*(3+tau)-18*(4-5*tau+tau^3)) ) / log(3/2) 

	# estimate of c (from approximated equation)	
	if (k < -1) {
		c <- -((k+(5/2)^(1+k)*k+k*tau+k*sqrt((2/5)^(-2*(1+k))+2^-k*5^(1+k)*(1-3*tau)+(1+tau)^2))/(4*tau))
	} else {
		c <- -(2/3)*(5/2-(2/5)^(-1-k))*(1+k-1/tau)+2/3*(-1+(5/2)^(1+k))*(-k-log(2)/log(1-tau))
	}
	
	# estimate of lambda (from exact exation, but a function of approximated k and c values)
	lambda <- L1*(-k)^(1+1/c)*gamma(1-1/k)/(gamma(1+1/c)*gamma(-1/c-1/k))	
	
	return(c(lambda=lambda, k=k, c=c))
}














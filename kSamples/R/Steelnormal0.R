Steelnormal0 <- function(sig0,sig,tau,U,ni,
		alternative=c("greater","two-sided")){
# this function computes the normal approximation of the p-value for the Steel test,
# based on the sizes ni = c(n1,...,nk) of the k treatment samples
# based on the observed Steel statistic U
# sig0, sig and tau are parameters required for the power evaluation. 
alternative <- match.arg(alternative)

k <- length(ni)
if(alternative=="greater"){ 
	funz <- function(z,k,sig0,sig,tau,S,ni){
		fac <- 1
		for(i in 1:k){
			fac <- fac * pnorm((S*tau[i]-ni[i]*sig0*z)/sig[i])
		}
		dnorm(z)*fac
	}
	N <- length(U)
	pval <- numeric(N)
	for(i in 1:N){
		pval[i] <- 1-integrate(funz,-Inf,Inf,k,sig0,sig,tau,U[i],ni)$value
	}
}


if(alternative=="two-sided"){ 

	funz <- function(z,k,sig0,sig,tau,S,ni){
			fac <- 1
			for(i in 1:k){
				fac <- fac * (pnorm((S*tau[i]-ni[i]*sig0*z)/sig[i])-
					pnorm((-S*tau[i]-ni[i]*sig0*z)/sig[i]))
			}
			dnorm(z)*fac
	}
	N <- length(U)
	pval <- numeric(N)
	for(i in 1:N){
		pval[i] <- 1-integrate(funz,-Inf,Inf,k,sig0,sig,tau,U[i],ni)$value
	}
}
pval
}

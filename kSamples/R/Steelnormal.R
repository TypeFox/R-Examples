Steelnormal <- function(mu,sig0,sig,tau,Wvec,ni,
		alternative=c("greater","less","two-sided"),continuity.corr=TRUE){
# this function computes the normal approximation of the p-value for the Steel test,
# based on the sizes ni = c(n1,...,nk) of the k treatment samples
# based on the vector of Mann-Whitney statistics comparing the treatment sample values (Y)
# against the common control sample values (X), Wvec consists of k such comparison statistics
# counting X_i < Y_j and 0.5 of X_i = Y_j.
# mu , sig0, sig, and, tau are parameters required for the power evaluation. 
alternative <- match.arg(alternative)
if(continuity.corr==TRUE){
	cont.corr <- .5
}else{
	cont.corr <- 0
}
k <- length(ni)
if(alternative=="greater"){ 
 Sx <- max((Wvec-mu)/tau)
 i0 <- min((1:k)[Sx == (Wvec-mu)/tau])
 S <- (Wvec[i0]-cont.corr-mu[i0])/tau[i0]
 funz <- function(z,k,sig0,sig,tau,S,ni){
		fac <- 1
		for(i in 1:k){
			fac <- fac * pnorm((S*tau[i]-ni[i]*sig0*z)/sig[i])
		}
		dnorm(z)*fac
	}
pval <- 1-integrate(funz,-Inf,Inf,k,sig0,sig,tau,S,ni)$value
}
if(alternative=="less"){ 
 Sx <- min((Wvec-mu)/tau)
 i0 <- min((1:k)[Sx == (Wvec-mu)/tau])
 S <- (Wvec[i0]+cont.corr-mu[i0])/tau[i0]
 funz <- function(z,k,sig0,sig,tau,S,ni){
		fac <- 1
		for(i in 1:k){
			fac <- fac * (1-pnorm((S*tau[i]-ni[i]*sig0*z)/sig[i]))
		}
		dnorm(z)*fac
	}
pval <- 1-integrate(funz,-Inf,Inf,k,sig0,sig,tau,S,ni)$value
}
if(alternative=="two-sided"){ 
 Sx <- max(abs(Wvec-mu)/tau)
 i0 <- min((1:k)[Sx == abs(Wvec-mu)/tau])
 S <- (abs(Wvec[i0]-mu[i0])-cont.corr)/tau[i0]
 funz <- function(z,k,sig0,sig,tau,S,ni){
		fac <- 1
		for(i in 1:k){
			fac <- fac * (pnorm((S*tau[i]-ni[i]*sig0*z)/sig[i])-
				pnorm((-S*tau[i]-ni[i]*sig0*z)/sig[i]))
		}
		dnorm(z)*fac
	}
pval <- 1-integrate(funz,-Inf,Inf,k,sig0,sig,tau,S,ni)$value
}
pval

}
# end of Steelnormal

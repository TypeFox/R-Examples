init.all <- function(nval,K,rand.start,mixture=FALSE) {
#
# Function init.all to create (rather arbitrary) initial values for
# tpm and Rho.
#

if(rand.start$tpm)
	tpm <- if(mixture) matrix(runif(K),K,K,byrow=TRUE)
			else matrix(runif(K*K),K,K)
else tpm <- matrix(1/K,K,K) + 1.5*diag(K)
tpm <- tpm/apply(tpm,1,sum)

if(rand.start$Rho) Rho <- matrix(runif(K*nval),K,nval)
else Rho <- matrix(1:(nval*K),K,nval)
Rho <- t(Rho/apply(Rho,1,sum))

list(tpm=tpm,Rho=Rho)
}

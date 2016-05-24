varma <- function(N, k=1, VAR = NULL, VMA = NULL, cov_matrix=diag(k), innov = NULL){
## Generates N observations of a k-vector ARMA process
##
## INPUT   N 		Number of time points
##	   k 		Dimension of the vectorial ARMA (optional, default is univariate)
##	   VAR		Array of VAR coefficient matrices (optional)
##	   VMA		Array of VMA coefficient matrices (optional)
##	   cov_matrix	Matrix of covariance between the innovations (optional, default is identity)
##         innov	Matrix of the innovations (optional, default is a gaussian process)
##         
## OUTPUT   		Vector containing the N observations of the k-vector ARMA process
##
##										Achard & Gannaz (2014)
##__________________________________________________________________________________________________________________

if((dim(cov_matrix)[1]!=k)|(dim(cov_matrix)[2]!=k)){ stop('The dimension of the covariance matrix mismatches') }

if(is.null(innov)==TRUE) innov<-matrix(rnorm(k*N),N)%*%chol(cov_matrix[1:k,1:k])

if(is.matrix(innov)){ 
	nn <- dim(innov)[1] 
	kk <- dim(innov)[2]
	}
else{ 
	nn <- length(innov) 
	kk <- 1
}
innov <- as.matrix(innov,dim=c(nn,kk))
if(kk!=k){ warning('The dimension of the innovations mismatches') }

if(is.null(VAR)==FALSE){ if(dim(VAR)[2]!=k){ stop('The dimension of AR does not match') }}
if(is.null(VMA)==FALSE){ if(dim(VMA)[2]!=k){ stop('The dimension of MA does not match') }}


# check if there are AR and MA components
nar <- 0
nma <- 0
if(is.null(VAR)==FALSE){
	nar <- dim(VAR)[3]
	if(is.na(nar)){ 
		nar <-1 
		VAR <- array(VAR,dim=c(k,k,1))
	}
}
if(is.null(VMA)==FALSE){
	nma <- dim(VMA)[3]
	if(is.na(nma)){ 
		nma <-1
		VMA <- array(VMA,dim=c(k,k,1))
	}
}


## simulate the data
u0 <- innov
if(nma>0){
	for(i in (nma+1):nn){
		for(lag in 1:nma){
			u0[i,] <- u0[i,] + innov[i-lag,]%*%VMA[,,lag]
		}
	}
}	
x <- u0
if(nar>0){
	for(i in (nar+1):nn){
		for(lag in 1:nar){
			x[i,] <- x[i,] - x[i-lag,]%*%VAR[,,lag]
		}
	}
}
x <- x[(nn-N+1):nn,]


return(x)
}


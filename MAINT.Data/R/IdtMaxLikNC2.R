GC2mLogLik <- function(t0,X,ue=NULL,maxsdratio=10000)  # Minus Log-likelihood for Gaussian model with Configuration 2
{
   # Wrapper function to be minimized in the ML estimation of a gaussian model assuming configuration 2 
   # (all correlations allowed except for mid-points and log-ranges between different interval variables) 
   #
   # Arguments:
   # 
   #  t0 - vector of parameters to be optimized. It should have as its components the non-null elements (row by row)
   #        of the lower-triangular Choleski decomposition of the covariance matrix among mid-points, followed by 
   #        the non-null elements for the log-ranges, and finally the non-null elements for the pairs mid-points
   #        log-ranges of the same variables 
   #
   #  Other arguments:
   #
   #  X  -  		Data matrix containing the mid-points in the first columns and the log-ranges in the following columns 
   #  ue -  		Matrix (with the same dimension as X) containing observation specific estimates of the X means 
   #	    		If NULL (default) the data is assumed to be centered
   #  maxsdratio -	Maxmimum allowale value for ratios of standard deviations before the covariance matrix is considered
   #                    to be numerically singular, and a (large) penalty is added to the negative log-likelihood 
   #
   #  Value: -1 times the log-likelihood 

	PenF <- 1E12					# large penalty for unfeasible parameters                             
	if ( any(!is.finite(t0)) ) return(PenF)   	# make sure that all parameters are valid real number

	p <- ncol(X)					# total number of variables (mid-points + log-ranges)
	q <- p/2					# number of interval variables

   	Sigmasr <- matrix(0.,nrow=p,ncol=p)
	intcomb <- q*(q+1)/2			# Number of possible combinations between pairs of interval variables
	for (i in 1:q) {
		for (j in 1:i) {		
			Sigmasr[i,j] <- t0[i*(i-1)/2+j]			# loadings among mid-points  
			Sigmasr[q+i,q+j] <- t0[intcomb+i*(i-1)/2+j]	# loadings among log-ranges
		}	 
		Sigmasr[q+i,i] <- t0[2*intcomb+i]	# loadings between mid-points and log-ranges of the same variables
	}
	for (i in 1:(q-1))		# non-free values of the loading matrix, required to ensure null-correlations
		for (j in (i+1):q)  {	# between mid-points and log-ranges of diferent interval variables	
			tempsum <- 0.		
			for (k in i:(j-1))  tempsum <- tempsum + Sigmasr[q+i,k]*Sigmasr[j,k]
			Sigmasr[q+i,j] <- -tempsum/Sigmasr[j,j]
		}
	Sigsd <- diag(Sigmasr)
	sdratio <- max(Sigsd)/min(Sigsd)
	if (sdratio > maxsdratio) return(PenF*(sdratio-maxsdratio))
	if (!is.null(ue)) X <- scale(X,center=ue,scale=FALSE) 

	-sum( apply(X, 1, ILogLikNC1, SigmaSrInv=forwardsolve(Sigmasr,diag(p)), const=-0.5*p*log(2*pi)-sum(log(diag(Sigmasr)))) )
}

GC2mLogLik.grad <- function(t0,X,ue=NULL)
{
    # gradient of GC2mLogLik 

	PenF <- 1E12					# large penalty for unfeasible parameters                             
	if ( any(!is.finite(t0)) ) return(PenF)   	# make sure that all parameters are valid real number

	p <- ncol(X)					# total number of variables (mid-points + log-ranges)
	n <- nrow(X)					# total number of observations
	q <- p/2					# number of interval variables

   	Sigmasr <- matrix(0.,nrow=p,ncol=p)
	intcomb <- q*(q+1)/2			# Number of possible combinations between pairs of interval variables
	difintcomb <- q*(q-1)/2			# Number of possible combinations between pairs of different interval variables
	for (i in 1:q) {
		for (j in 1:i) {		
			Sigmasr[i,j] <- t0[i*(i-1)/2+j]		# loadings among mid-points  
			Sigmasr[q+i,q+j] <- t0[intcomb+i*(i-1)/2+j]	# loadings among log-ranges
		}	 
		Sigmasr[q+i,i] <- t0[2*intcomb+i]	# loadings between mid-points and log-ranges of the same variables
	}
	for (i in 1:(q-1))		# non-free values of the loading matrix, required to ensure null-correlations
		for (j in (i+1):q)  {	# between mid-points and log-ranges of diferent interval variables	
			tempsum <- 0.		
			for (k in i:(j-1))  tempsum <- tempsum + Sigmasr[q+i,k]*Sigmasr[j,k]
			Sigmasr[q+i,j] <- -tempsum/Sigmasr[j,j]
		}
        if (is.null(ue)) {
		if (!is.matrix(X)) X <- as.matrix(X)
		V = t(X) %*% X / n
	}
	else  {
		Xdev = as.matrix(X-ue)
        	V = t(Xdev) %*% Xdev / n
	}
	SigmaInv <- chol2inv(t(Sigmasr))
	fixpgrad = matrix(0.,nrow=2*q,ncol=2*q)
	for (a in 1:(q-1))  for (b in q:(a+1)) {
		fixpgrad[q+a,b] = fixpgrad[q+a,b] + fL.grad(dfx=Sigmagrad,L=Sigmasr,l1=q+a,l2=b,totald=TRUE,L2=SigmaInv,n=n,V=V) 
		if (b > a+1) for (k in (a+1):(b-1)) fixpgrad[q+a,k] = fixpgrad[q+a,k] - Sigmasr[b,k]/Sigmasr[b,b] * fixpgrad[q+a,b]
	}

 	grad = array(dim=2*intcomb+q)
	for (i in 1:q)  { 
		for (j in 1:i) {		
			grad[i*(i-1)/2+j] = fL.grad(dfx=Sigmagrad,L=Sigmasr,l1=i,l2=j,totald=TRUE,L2=SigmaInv,n=n,V=V) 
			if (j<i) for (a in 1:(i-1)) grad[i*(i-1)/2+j] = grad[i*(i-1)/2+j] - Sigmasr[q+a,j]/Sigmasr[i,i] * fixpgrad[q+a,i]
			else if (j>1) for (a in 1:(j-1)) for (k in a:(j-1))
				grad[j*(j-1)/2+j] = grad[j*(j-1)/2+j] + Sigmasr[q+a,k]*Sigmasr[j,k]/Sigmasr[j,j]^2 * fixpgrad[q+a,j]
			grad[intcomb+i*(i-1)/2+j] = fL.grad(dfx=Sigmagrad,L=Sigmasr,l1=q+i,l2=q+j,totald=TRUE,L2=SigmaInv,n=n,V=V) 
		}
		grad[2*intcomb+i] = fL.grad(dfx=Sigmagrad,L=Sigmasr,l1=q+i,l2=i,totald=TRUE,L2=SigmaInv,n=n,V=V)
		if (i<q) for (b in (i+1):q) grad[2*intcomb+i] = grad[2*intcomb+i] - Sigmasr[b,i]/Sigmasr[b,b] * fixpgrad[q+i,b]
	}
	-matrix(grad,1,length(grad))
}

GC2LogLikC.grad <- function(t0,X,ue=NULL)   
# Gradient of observation contributions for Log-likelihood of Gaussian model with Configuration 2
{
	PenF <- 1E12					# large penalty for unfeasible parameters                             
	if ( any(!is.finite(t0)) ) return(PenF)   	# make sure that all parameters are valid real number

	p <- ncol(X)					# total number of variables (mid-points + log-ranges)
	n <- nrow(X)					# total number of observations
	q <- p/2					# number of interval variables

   	Sigmasr <- matrix(0.,nrow=p,ncol=p)
	intcomb <- q*(q+1)/2			# Number of possible combinations between pairs of interval variables
	difintcomb <- q*(q-1)/2			# Number of possible combinations between pairs of different interval variables
	for (i in 1:q) {
		for (j in 1:i) {		
			Sigmasr[i,j] <- t0[i*(i-1)/2+j]		# loadings among mid-points  
			Sigmasr[q+i,q+j] <- t0[intcomb+i*(i-1)/2+j]	# loadings among log-ranges
		}	 
		Sigmasr[q+i,i] <- t0[2*intcomb+i]	# loadings between mid-points and log-ranges of the same variables
	}
	for (i in 1:(q-1))		# non-free values of the loading matrix, required to ensure null-correlations
		for (j in (i+1):q)  {	# between mid-points and log-ranges of diferent interval variables	
			tempsum <- 0.		
			for (k in i:(j-1))  tempsum <- tempsum + Sigmasr[q+i,k]*Sigmasr[j,k]
			Sigmasr[q+i,j] <- -tempsum/Sigmasr[j,j]
		}
        if (is.null(ue)) { if (!is.matrix(X)) X <- as.matrix(X) }
	else  X <- as.matrix(X-ue)
        Vcnt <- apply(X,1,function(v) outer(v,v)/n)
	dim(Vcnt) <- c(p,p,n)
	SigmaInv <- chol2inv(t(Sigmasr))
	fixpgradc = array(0.,dim=c(p,p,n))
	for (a in 1:(q-1))  for (b in q:(a+1)) {
		fixpgradc[q+a,b,] = fixpgradc[q+a,b,] + apply(Vcnt,3,SigmasrgradCnt,L=Sigmasr,l1=q+a,l2=b,L2=SigmaInv,n=n) 
		if (b > a+1) for (k in (a+1):(b-1)) fixpgradc[q+a,k,] = fixpgradc[q+a,k,] - Sigmasr[b,k]/Sigmasr[b,b] * fixpgradc[q+a,b,]
	}
 	gradc = matrix(nrow=2*intcomb+q,ncol=n)
	for (i in 1:q)  { 
		for (j in 1:i) {		
			gradc[i*(i-1)/2+j,] =  apply(Vcnt,3,SigmasrgradCnt,L=Sigmasr,l1=i,l2=j,L2=SigmaInv,n=n) 
			if (j<i) for (a in 1:(i-1)) gradc[i*(i-1)/2+j,] = gradc[i*(i-1)/2+j,] - Sigmasr[q+a,j]/Sigmasr[i,i] * fixpgradc[q+a,i,]
			else if (j>1) for (a in 1:(j-1)) for (k in a:(j-1))
				gradc[j*(j-1)/2+j,] = gradc[j*(j-1)/2+j,] + Sigmasr[q+a,k]*Sigmasr[j,k]/Sigmasr[j,j]^2 * fixpgradc[q+a,j,]
			gradc[intcomb+i*(i-1)/2+j,] =  apply(Vcnt,3,SigmasrgradCnt,L=Sigmasr,l1=q+i,l2=q+j,L2=SigmaInv,n=n) 
		}
		gradc[2*intcomb+i,] = apply(Vcnt,3,SigmasrgradCnt,L=Sigmasr,l1=q+i,l2=i,L2=SigmaInv,n=n) 
		if (i<q) for (b in (i+1):q) gradc[2*intcomb+i,] = gradc[2*intcomb+i,] - Sigmasr[b,i]/Sigmasr[b,b] * fixpgradc[q+i,b,]
	}
	t(gradc)
}

initparconf2 <- function(Data,n,q)	
{					 
# 	Find initial parameter estimates in configuration 2 (all correlations allowed except for mid-points and log-ranges 
#	between different interval variables) by replacing values in the Choleski decomposition of the within covariance 
#	matrix in order to inforce the required null-correlations

	EPS <- 1E-9
	intcomb <- q*(q+1)/2			# Number of possible combinations between pairs of interval variables
	initpar <- array(dim=2*intcomb+q)
	S <- var(Data)*(n-1)/n	
	Segv <- eigen(S,symmetric=TRUE,only.values=TRUE)$values	
	if (Segv[2*q]<EPS) Sigmasr <- diag(1,2*q)			
	else Sigmasr <- t(chol(S))	
	cnt <- 1	
	for (i in 1:q) {
		for (j in 1:i)  {
			initpar[cnt] <- Sigmasr[i,j] 
			initpar[intcomb+cnt] <- Sigmasr[q+i,q+j]
			cnt <- cnt + 1
		}
		initpar[2*intcomb+i] <- Sigmasr[q+i,i]
	}
	initpar
} 

C2GetCov <- function(par,OStdv,q,tol=1E-8)

#  Gets the covariance matrix from the optimal Configuration 2 optimization applied to scaled data

#  Arguments

#  par    - The optimal parameter values from the (scaled data) Configuration 2 optimization
#  OStdv  - Matrix with the outer product of the standard deviations
#  q      - Number of integer variables
  
{
	L <- matrix(0.,2*q,2*q)
	intcomb <- q*(q+1)/2		# Number of possible combinations between pairs of interval variables
	cnt <- 1	
	for (i in 1:q) {
		for (j in 1:i)  {
			L[i,j] <- par[cnt]  
			L[q+i,q+j] <- par[intcomb+cnt]  
			cnt <- cnt + 1
		}
		L[q+i,i] <- par[2*intcomb+i]  
	}
	for (i in 1:q)  {					# non-free values of the loading matrix, required to 
		if (i>1) for (j in 1:(i-1)) L[q+i,j] <- 0. 	# ensure null-correlations between mid-points and 	
		if (i<q) for (j in (i+1):q)  {			# log-ranges of diferent interval variables	
			tempsum <- 0.		
			for (k in i:(j-1))  tempsum <- tempsum + L[q+i,k]*L[j,k]
			L[q+i,j] <- -tempsum/L[j,j]
		}
	}
	StdSigma <- L %*% t(L)
	StdSigma[abs(StdSigma)<tol]  <- 0.
	OStdv * StdSigma
}

Cnf2MaxLik <- function(Data,sd0=0.1,maxrepet=3,maxnoimprov=20,maxreplic=200,iter=10000,eval=100000,EPS=1E-6)
{
   #  Note  -  The Data argument should be a matrix containing the mid-points in the first columns and the log-ranges in the following columns 

	n <- nrow(Data)			    # number of observations
	p <- ncol(Data)			    # total number of variables (mid-points + log-ranges)
	q <- p/2			    # number of interval variables
	intcomb <- q*(q+1)/2	   	    # Number of possible combinations between pairs of interval variables
	npar <- 2*intcomb + q	   	    # Total number of parameters to optimize
	initpar <- initparconf2(Data,n,q)	
	parsd <- rep(sd0,npar)      	    # standard deviation hyper-parameters - used to generate random starting points
	parlb <- rep(-Inf,npar)    	    # vector of lower bounds
	for (i in 1:q)  { 
	  parlb[i*(i-1)/2+i] <- EPS 	    # diagonal elements of Choleski decomposition must be positive
	  parlb[intcomb +i*(i-1)/2+i] <- EPS
	}
	res <- RepLOptim(initpar,parsd,fr=GC2mLogLik,gr=GC2mLogLik.grad,
			maxrepet=maxrepet,maxnoimprov=maxnoimprov,maxreplic=maxreplic,
			method="nlminb",maxiter=iter,maxeval=eval,lower=parlb,X=Data)
	list(lnLik=-res$val,SigmaSr=res$par,optres=res)
}

Sigmagrad <- function(j1,j2,n,V,L2)
{
   grad <- n * ( sum(outer(L2[,j1],L2[,j2])*V) - L2[j1,j2] )
   if (j1==j2)  grad <- grad / 2
   grad #  return(grad) 
}

SigmasrgradCnt <- function(VCnt,L2,L,l1,l2,n)  fL.grad(dfx=SigmagradCnt,L=L,l1=l1,l2=l2,totald=TRUE,L2=L2,n=n,VCnt=VCnt)

SigmagradCnt <- function(VCnt,InvMat,j1,j2,n)
{
    grad <- n*sum(outer(InvMat[,j1],InvMat[,j2])*VCnt) - InvMat[j1,j2] 
    if (j1==j2)  grad <- grad / 2
    grad #  return(grad) 
}


VCovGC2LogLikC.grad <- function(t0,X,ue=NULL)   
# Variance-Covariance Gradient of observation contributions for Log-likelihood of Gaussian model with Configuration 2
{
	PenF <- 1E12					# large penalty for unfeasible parameters                             
	if ( any(!is.finite(t0)) ) return(PenF)   	# make sure that all parameters are valid real number

	p <- ncol(X)					# total number of variables (mid-points + log-ranges)
	n <- nrow(X)					# total number of observations
	q <- p/2					# number of interval variables

   	Sigma <- matrix(0.,nrow=p,ncol=p)
	intcomb <- q*(q+1)/2			# Number of possible combinations between pairs of interval variables
	difintcomb <- q*(q-1)/2			# Number of possible combinations between pairs of different interval variables
	for (i in 1:q) {
		for (j in 1:i) {		
			Sigma[i,j] <- Sigma[j,i] <- t0[i*(i-1)/2+j]			# covariances among mid-points  
			Sigma[q+i,q+j] <- Sigma[q+j,q+i] <- t0[intcomb+i*(i-1)/2+j]	# covariances among log-ranges
		}	 
		Sigma[q+i,i] <- Sigma[i,q+i] <- t0[2*intcomb+i]		# covariances between mid-points and log-ranges of the same variables
	}
        if (is.null(ue)) { if (!is.matrix(X)) X <- as.matrix(X) }
	else  X <- as.matrix(X-ue)
        Vcnt <- apply(X,1,function(v) outer(v,v)/n)
	dim(Vcnt) <- c(p,p,n)
	SigmaInv <- solve(Sigma)
 	gradc = matrix(nrow=2*intcomb+q,ncol=n)
	for (i in 1:q)  { 
		for (j in 1:i) {		
			gradc[i*(i-1)/2+j,] =  apply(Vcnt,3,SigmagradCnt,j1=i,j2=j,InvMat=SigmaInv,n=n) 
			gradc[intcomb+i*(i-1)/2+j,] =  apply(Vcnt,3,SigmagradCnt,j1=q+i,j2=q+j,InvMat=SigmaInv,n=n) 
		}
		gradc[2*intcomb+i,] = apply(Vcnt,3,SigmagradCnt,j1=q+i,j2=i,InvMat=SigmaInv,n=n) 
	}
	t(gradc)
}

CovtoParC2 <- function(Sigma,q,Sqrt)	
{					 
# 	Gets the vectorized form of covariance matrix according to Configuration 2

	intcomb <- q*(q+1)/2	   	    # Number of possible combinations between pairs of interval variables
	par <- array(dim=2*intcomb+q)
	if (Sqrt==TRUE) mat <- t(chol(Sigma))
	else mat <- Sigma
	cnt <- 1	
	for (i in 1:q) {
		for (j in 1:i)  {
			par[cnt] <- mat[i,j] 
			par[intcomb+cnt] <- mat[q+i,q+j]
			cnt <- cnt + 1
		}
		par[2*intcomb+i] <- mat[q+i,i]
	}
	par
} 

PartoMatC2 <- function(par,q,initval)

#  Converts a vectorized form of parameters (or standard errors) according to Configuration 2, to ist matrix form


{
	mat <- matrix(initval,2*q,2*q)
	intcomb <- q*(q+1)/2		# Number of possible combinations between pairs of interval variables
	cnt <- 1	
	for (i in 1:q) {
		for (j in 1:i)  {
			mat[i,j] <- mat[j,i] <- par[cnt]  
			mat[q+i,q+j] <- mat[q+j,q+i] <- par[intcomb+cnt]  
			cnt <- cnt + 1
		}
		mat[q+i,i] <- mat[i,q+i] <- par[2*intcomb+i]  
	}
	mat
}

C2GetCovStderr <- function(mleSigmaC2,Data,q,ue=NULL)

#  Arguments

#  mleSigmaC2    - The maximum likelihood estimators of the covariance matrix under configuration C2
#  Data      	 - Data frame with the MidPoints followed by the Log-Ranges of the Interval Data set.
#                  When ue is set to NULL (default) Data is assumed to be centered.
#  q             - Number of Integer Varaibles
#  ue		 - Vector with the MidPoints and Log-Ranges means (or NULL if Data is centered)
  
{
	p <- 2*q
	n <- nrow(Data)
	if (!is.null(ue)) Data <- scale(Data,center=ue,scale=FALSE)
	par <- CovtoParC2(mleSigmaC2,q,Sqrt=FALSE)
	npar <- length(par)	
	gradc <- VCovGC2LogLikC.grad(par,Data)
	OtProd <- apply(gradc,1,function(v) outer(v,v))
	dim(OtProd) <- c(npar,npar,n)
	HessAp <- matrix(apply(OtProd,c(1,2),sum),nrow=npar,ncol=npar,byrow=FALSE)
	stderr <- sqrt(diag(solve(HessAp)))
	PartoMatC2(stderr,q,NA)
}


###################################
# Fonctions called by PTCMestimBF #
###################################

# obtainLambda
##############
# Solve Equations (4) and (6) in Ma and Yin (2008) for lambda
# Inputs :
  # - lambda: the Lagrange multiplier
  # - betaVec: vector (of length R+P+1) containing the regression coefficients
  # - m: the number of distinct failure times
  # - R: the number of covariates without measurement error
  # - P: the number of mismeasured covariates
  # - BOrd: vector (of length m) containing the ordered distinct failure times
  # - ObsTime: vector (of length n) containing all the ordered observed times
  # - Xc: matrix (of dimensions n x R) containing the covariates without measurement error, ordered according to ObsTime
  # - V: matrix (of dimensions n x P) containing the mismeasured covariates, ordered according to ObsTime
  # - OmegaMatEstim: variance-covariance matrix (of dimensions (R+P+1)x(R+P+1)) of the measurement error
  # - n: number of observations
  
obtainLambda <- function(lambda,betaVec,m,R,P,BOrd,ObsTime,Xc,V,OmegaMatEstim,n) 
{
	W <- cbind(1,V,Xc)
	bvb <- 0.5*matrix(nrow=1,ncol=R+P+1,data=betaVec)%*%OmegaMatEstim%*%betaVec
	expb <- exp(W%*%betaVec-c(bvb))
	useful <- matrix(nrow=m, ncol=n, byrow=TRUE, data= (rep(BOrd,rep(n,m))<=rep(ObsTime[,1],m) ) )
	useful2 <- matrix(nrow=m,ncol=n,byrow=TRUE,data=(ObsTime[,2]<2))
	useful3 <- (useful*useful2)%*%expb+(dim(W)[1])*lambda
	
	sum(1/useful3)-1
}

# obtainPs
#########
# Solve Equation (6) in Ma and Yin (2008) for the p_(i)'s
# Inputs:
  # - lambda: the Lagrande multiplier
  # - betaVec: vector (of length R+P+1) containing the regression coefficients
  # - m: the number of distinct failure times
  # - R: the number of covariates without measurement error
  # - P: the number of mismeasured covariates
  # - BOrd: vector (of length m) containing the ordered distinct failure times
  # - ObsTime: vector (of length n) containing all the ordered observed times
  # - Xc: matrix (of dimensions n x R) containing the covariates without measurement error, ordered according to ObsTime
  # - V: matrix (of dimensions n x P) containing the mismeasured covariates, ordered according to ObsTime
  # - OmegaMatEstim: variance-covariance matrix (of dimensions (R+P+1)x(R+P+1)) of the measurement error
  # - n: number of observations
  
obtainPs <- function(lambda,betaVec,m,R,P,BOrd,ObsTime,Xc,V,OmegaMatEstim,n)
{
	W <- cbind(1,V,Xc)
	bvb <- 0.5*matrix(nrow=1,ncol=R+P+1,data=betaVec)%*%OmegaMatEstim%*%betaVec
	expb <- exp(W%*%betaVec-c(bvb))
	useful <- matrix(nrow=m, ncol=n, byrow=TRUE, data= (rep(BOrd,rep(n,m))<=rep(ObsTime[,1],m) ) )
	useful2 <- matrix(nrow=m,ncol=n,byrow=TRUE,data=(ObsTime[,2]<2))	
	useful3 <- (useful*useful2)%*%expb+(dim(W)[1])*lambda
	
	1/useful3
}

# obtainBeta
############
# Solve Equation (7) in Ma and Yin (2008) for beta
# Inputs:
  # - betaVec: vector (of length R+P+1) containing the regression coefficients
  # - pOrd: vector (of length n) containing the steps of the baseline cdf
  # - m: the number of distinct failure times
  # - R: the number of covariates without measurement error
  # - P: the number of mismeasured covariates
  # - BOrd: vector (of length m) containing the ordered distinct failure times
  # - ObsTime: matrix (of dimensions n x 2) containing all the ordered observed times in the 1st column, and the status of each observation in the 2nd column (0=censored, 1=failure, 2=cured)
  # - Xc: matrix (of dimensions n x R) containing the covariates without measurement error, ordered according to ObsTime
  # - V: matrix (of dimensions n x P) containing the mismeasured covariates, ordered according to ObsTime
  # - OmegaMatEstim: variance-covariance matrix (of dimensions (R+P+1)x(R+P+1)) of the measurement error

obtainBeta <- function(betaVec,pOrd,m,R,P,BOrd,ObsTime,Xc,V,OmegaMatEstim)
{
	W <- cbind(1,V,Xc)
	n <- dim(ObsTime)[1]
	bvb <- 0.5*matrix(nrow=1,ncol=R+P+1,data=betaVec)%*%OmegaMatEstim%*%betaVec
	expb <- exp(W%*%betaVec-c(bvb))

	useful <- matrix(nrow=m, ncol=n, byrow=TRUE, data= (rep(BOrd,rep(n,m))<=rep(ObsTime[,1],m) ) )
	useful1 <- (ObsTime[,2]>0.5 & ObsTime[,2]<1.5)
	useful2 <- t(useful)%*%pOrd
	useful3bis <- matrix(nrow=n,ncol=R+P+1,data=useful1,byrow=FALSE)
	useful3ter <- matrix(nrow=n,ncol=R+P+1,data=useful2*expb,byrow=FALSE)
	colSums(useful3bis*W-useful3ter*(W-matrix(nrow=n,ncol=R+P+1,data=OmegaMatEstim%*%betaVec,byrow=TRUE)))
}

# MaYinEstimation
#################
# Solve Equations (4), (6) and (7) of Ma and Yin (2008) by a backfitting approach
# Inputs:
  # - mm: the number of distinct failure times
  # - R: the number of covariates without measurement error
  # - P: the number of mismeasured covariates
  # - BBOrd: vector (of length m) containing the ordered distinct failure times
  # - ObsTime: matrix (of dimensions n x 2) containing all the ordered observed times in the 1st column, and the status of each observation in the 2nd column (0=censored, 1=failure, 2=cured)
  # - Xc: matrix (of dimensions n x R) containing the covariates without measurement error, ordered according to ObsTime
  # - V: matrix (of dimensions n x P) containing the mismeasured covariates, ordered according to ObsTime
  # - OmegaMatEstim: variance-covariance matrix (of dimensions (R+P+1)x(R+P+1)) of the measurement error  
  # - init: vector (of length R+P+1) containing initial values for the regression coefficients beta
  # - nBack: maximum number of iterations in the backfitting approach
  # - eps: convergence criterion
  # - multTimeMYestim: when multiplied by 300, gives the maximal number of seconds during which function nleqslv is allowed to run before being stopped.
  # - n: number of observations

MaYinEstimation <- function(mm,R,P,BBOrd,OObsTime,Xc,V,OmegaMatEstim,init,nBack,eps,multTimeMYestim,n)
{
	newLambda <- rep(NA,nBack)
	newpOrd <- matrix(nrow=nBack,ncol=mm)
	newbetaVec <- matrix(nrow=nBack,ncol=1+R+P)
	
	newLambda[1] <- 0.05#runif(1)
	newpOrd[1,] <- 1/mm
	newbetaVec[1,] <- init
	
	endK <- 0
	flag <- 0
	flag1 <- rep(0,nBack)
	flag2 <- rep(0,nBack)	
	
	for(k in 2:nBack)
	{
		comparison <- c(newLambda[k-1],newpOrd[k-1,],newbetaVec[k-1,])
		
		time1 <- proc.time()[[3]]
	
		j=0
		# Step 1a : solve for lambda
		while(flag1[k]<0.5) # Function nleqslv often requires multiple runs (with different starting values) before achieving convergence
		{
			j <- j+1
			if((proc.time()[[3]]-time1)<300*multTimeMYestim)
			{
				lambdastart <- (j==1)*newLambda[k-1]+(j!=1)*runif(1)  
				getLambda <- nleqslv(x=lambdastart, method="Broyden",fn=obtainLambda,betaVec=newbetaVec[k-1,],m=mm,R=R,P=P,BOrd=BBOrd,ObsTime=OObsTime,Xc=Xc,V=V,OmegaMatEstim=OmegaMatEstim,n=n)
				if (getLambda$termcd==1) # If convergence ok
				{
					newLambda[k] <- getLambda$x 
					flag1[k] <- 1
				} else
				{
					flag1[k] <- 0
				}
			}else{flag1[k]=2}  # flag1=1 if ok, 2 if convergence not reached after some predefined time.
		}		
		
		# Step 1b : solve for the p_i's
		newpOrd[k,] <- obtainPs(lambda=newLambda[k],betaVec=newbetaVec[k-1,],m=mm,R=R,P=P,BOrd=BBOrd,ObsTime=OObsTime,Xc=Xc,V=V,OmegaMatEstim=OmegaMatEstim,n=n)
	
		# Step 3 : solve for beta
		j=0
		while(flag2[k]<0.5) # Function nleqslv often requires multiple runs (with different starting values) before achieving convergence
		{
			j=j+1
			if((proc.time()[[3]]-time1)<300*multTimeMYestim)
			{
				betastart <- (j==1)*newbetaVec[k-1,]+(j!=1)*0.1*rnorm(R+1+P) 
				getBeta <- nleqslv(x=betastart, fn=obtainBeta,control=list(maxit=1000,xtol=1e-6),pOrd=newpOrd[k,],m=mm,R=R,P=P,BOrd=BBOrd,ObsTime=OObsTime,Xc=Xc,V=V,OmegaMatEstim=OmegaMatEstim)
				                                       
				if (getBeta$termcd==1) # If convergence ok
				{
					newbetaVec[k,] <- getBeta$x 
					flag2[k] <- 1
				} else
				{
					flag2[k] <- 0
				}
			}else{flag2[k]=2} # flag2=1 if ok, 2 if convergence not reached after some predefined time.
		}
		fct1 <- getLambda$fvec
		fct2 <- getBeta$fvec
		crit1 <- sqrt(crossprod(comparison-c(newLambda[k],newpOrd[k,],newbetaVec[k,])))
		if(crit1<eps)
		{
			crit2 <- sqrt(fct1^2)
			crit3 <- sqrt(crossprod(fct2))
			if(crit2<eps & crit3<eps) 
			{
				endK <- k; flag <- 1 ; break
			}
		}
	}
	if(k==nBack)
	{
		cat("nBack(",nBack,") reached","\n")	
		flag <- 2
		endK <- k
	}
	if(flag==2)  stop("NON CONVERGENCE")
	else return(list(estimMYI=newbetaVec[endK,],mp=newpOrd[endK,],flag=flag,flag1=any(flag1==2),flag2=any(flag2==2),endK=endK))
}


# computeVariance
#################

gett1 <- function(x,bet,sigma)
{
	n <- dim(x)[1]
	t1 <- exp(x%*%bet-c(t(bet)%*%sigma%*%bet/2))
}

gett2 <- function(x,bet,sigma,nbParam,n)
{
	t2 <- matrix(nrow=1+nbParam,ncol=n)
	t2 <- t(x)-c(sigma%*%bet)
	t2
}

gett3 <- function(x)
{
	t3 <- t(x)
	t3
}

getb1 <- function(d,y,t1,t2,n,nbParam)
{
	t(((matrix(nrow=n,ncol=n,byrow=TRUE,data=rep((d==2),n))+matrix(nrow=n,ncol=n,byrow=TRUE,data=(rep(y,rep(n,n))<=rep(y,n))))>0)%*%(matrix(nrow=n,ncol=1+nbParam,byrow=FALSE,data=t1)*t(t2)))
}

getb2 <- function(d,fy,b1,n) 
{
	(d[1]==1)*fy[1]*b1[,1]+b1[,2:n]%*%((fy[2:n]-fy[1:(n-1)])*(d[2:n]==1))
}

getc1 <- function(t1,fy,d)
{
	c(t1)%*%(c(fy)^(d<2))-d%*%(d<2)
}

getc2 <- function(t1,y,d,n)
{
	((matrix(nrow=n,ncol=n,byrow=TRUE,data=(rep(y,rep(n,n))<=rep(y,n)))+matrix(nrow=n,ncol=n,byrow=TRUE,data=rep((d==2),n)))>0)%*%t1
}

getc12 <- function(c1,c2)
{
	c12 <- 1/(c(c1)-c2)
}

getb3 <- function(d,c12,fy,b1,b2,n,nbParam) 
{
	tmpd <- (d[1]==1)*c12[1]*fy[1]+((d[2:n]==1)*c12[2:n])%*%(fy[2:n]-fy[1:(n-1)])
	tmpu21 <- (d[1]==1)*c12[1]*fy[1]*(b1[,1]-b2)+t( c((d[2:n]==1)*c12[1]*(fy[2:n]-fy[1:(n-1)]))%*%(t(b1[,2:n])-matrix(nrow=n-1,ncol=1+nbParam,byrow=TRUE,data=rep(b2,n-1))) )
	tmpu22 <- (d[1]==1)*c12[1]*fy[1]*(b1[,1]-b2)+t( c((d[2:n]==1)*c12[2:n]*(fy[2:n]-fy[1:(n-1)]))%*%(t(b1[,2:n])-matrix(nrow=n-1,ncol=1+nbParam,byrow=TRUE,data=rep(b2,n-1))) )
	c(tmpu22)/tmpd
}

getb4 <- function(c12,b1,b2,b3,nbParam,n) 
{
	matrix(nrow=1+nbParam,ncol=n,byrow=TRUE,data=c12)*(b1-matrix(nrow=1+nbParam,ncol=n,byrow=FALSE,data=b2)-matrix(nrow=1+nbParam,ncol=n,byrow=FALSE,data=b3))
}

geta1 <- function(sigma,t1,fy,x,bet,d,nbParam,n)
{
	a1 <- matrix(nrow=1+nbParam,ncol=1+nbParam,0)
	tmp <- matrix(nrow=1+nbParam,ncol=1+nbParam,0)
	
	for (i in 1:n)
	{
		tmp <- sigma-matrix(nrow=1+nbParam,ncol=1,data=x[i,]-sigma%*%bet)%*%matrix(ncol=1+nbParam,nrow=1,data=x[i,]-sigma%*%bet)
	
		a1 <- a1+fy[i]^(d[i]!=2)*t1[i]*tmp
	}
	a1
}

geta2 <- function(t1,t2,d,b4,fy,nbParam,n)
{
	a2 <- matrix(nrow=1+nbParam,ncol=1+nbParam,data=0)
	tmp1 <- rep(0,1+nbParam)
	tmp2 <- rep(0,1+nbParam)
	
	for (i in 1:n)
	{
		tmp1 <- t1[i]*t2[,i]
		tmp2 <- rep(0,1+nbParam)
		tmp2 <- (d[1]==1)*fy[1]*b4[,1]+(i>1)* b4[,2:n]%*%((d[2:n]==1)*(fy[2:n]-fy[1:(n-1)])*(c(2:n)<=i))
		a2 <- a2+tmp1%*%matrix(nrow=1,ncol=1+nbParam,data=tmp2)
	}
	a2
}

getsb <- function(t1,t2,fy,d,nbParam,n) 
{
	-t2*matrix(nrow=1+nbParam,ncol=n,byrow=TRUE, data=t1)*(matrix(nrow=1+nbParam,ncol=n,byrow=TRUE, data=fy))^matrix(nrow=1+nbParam,ncol=n,byrow=TRUE,data=(d!=2))+t2*matrix(nrow=1+nbParam,ncol=n,byrow=TRUE,data=(d!=2))*matrix(nrow=1+nbParam,ncol=n,byrow=TRUE,data=d)

}

getsf <- function(d,b4,fy,t1,nbParam,n)
{
	sftmp <- matrix(nrow=1+nbParam,ncol=n,byrow=FALSE, data=b4[,1]*fy[1]*(d[1]==1))+(b4[,2:n]*matrix(nrow=1+nbParam,ncol=n-1,byrow=TRUE,data=(d[2:n]==1))*matrix(nrow=1+nbParam,ncol=n-1,byrow=TRUE,data=fy[2:n]-fy[1:(n-1)]))%*%matrix(nrow=n-1,ncol=n,byrow=TRUE,data=(rep(2:n,rep(n,n-1))<=rep(1:n,n-1)))
	-sftmp*matrix(nrow=1+nbParam,ncol=n,byrow=TRUE,data=t1)+matrix(nrow=1+nbParam,ncol=n,byrow=TRUE,data=(d!=2))*matrix(nrow=1+nbParam,ncol=n,byrow=TRUE,data=d)*b4

}

getVar <- function(a1,a2,sb,sf)
{
	A <- t(a1-a2)
	B <- (sb+sf)%*%t(sb+sf)
	Var <- solve(A)%*%B%*%t(solve(A))
}


# computeVariance
#################
# Computes the asymptotic variance given in Theorem 3 in Ma and Yin (2008)
# Inputs:
  # - x: matrix (of dimensions n x (R+P+1)) containing a column of 1s, and the covariates 
  # - bet: vector (of length R+P+1) containing the regression coefficients
  # - sigma: variance-covariance matrix (of dimensions (R+P+1) x (R+P+1)) of the measurement errors
  # - d: vector (of length n) containing the status of each observation (0=censored, 1=failure, 2=cured)
  # - y: vector (of length n) containing the ordered observed times
  # - fy: vector (of length n) the estimated baseline cdf (ordered)
  # - n: number of observations
  # - nbParam: total number of regression parameters
  
computeVariance <- function(x,bet,sigma,d,y,fy,n,nbParam)
{
	t1 <- gett1(x,bet,sigma)
	t2 <- gett2(x,bet,sigma,nbParam,n)
	t3 <- gett3(x)
	b1 <- getb1(d,y,t1,t2,n,nbParam)
	b2 <- getb2(d,fy,b1,n) 
	c1 <- getc1(t1,fy,d)
	c2 <- getc2(t1,y,d,n)
	c12 <- getc12(c1,c2)
	b3 <- getb3(d,c12,fy,b1,b2,n,nbParam)
	b4 <- getb4(c12,b1,b2,b3,nbParam,n)
	a1 <- geta1(sigma,t1,fy,x,bet,d,nbParam,n) 
	a2 <- geta2(t1,t2,d,b4,fy,nbParam,n) 
	sb <- getsb(t1,t2,fy,d,nbParam,n)
	sf <- getsf(d,b4,fy,t1,nbParam,n) 
	VarMY <- getVar(a1,a2,sb,sf)
}	


# dataPreparation
#################
# Prepares the data before using a correction method: assign a cured status to some individuals, order the data according to the observed times
# Inputs:
  # - Dat: matrix (of dimensions n x (R+P+2)) containing, in the first 2 columns, the observed times and status, and then the covariates (with no column of 1's)
  # - P: number of mismeasured covariates
  # - R: number of covariates without measurement error 
	
dataPreparation <- function(Dat,P,R)
{
	OObsTime <- Dat[,1:2]
	V <- array(dim=c(dim(Dat)[1],P),Dat[,3:(2+P)])
	Xc <- array(dim=c(dim(Dat)[1],R),Dat[,(P+3):(P+2+R)])

	Tmax <- max(OObsTime[OObsTime[,2]==1,1])
	d <- OObsTime[,2]
	d[OObsTime[,2]==0 & OObsTime[,1]>Tmax] <- 2

	indicI <- rep(NA,3)
	indicI[1] <- sum(d==0)
	indicI[2] <- sum(d==1)
	indicI[3] <- sum(d==2) # The individuals considered as cured for the estimation procedure
 
	mm <- dim(table(OObsTime[(OObsTime[,2]==1),1])) # Number of distinct failure times (not taking into account the infinite failure times)	
	BBOrd <- sort(as.numeric(c(names(table(OObsTime[(OObsTime[,2]==1),1]))))) # ordered B_i, for B_i corresponding to event times (not censoring times) : B_(1)<B_(2)<...<B_(m)

	orderTime <- order(OObsTime[,1])
	OObsTime[,2] <- OObsTime[orderTime,2]
	V <- V[orderTime,]
	Xc <- Xc[orderTime,]
	OObsTime[,1] <- OObsTime[orderTime,1]

	return(list(OObsTime=OObsTime,V=V,Xc=Xc,mm=mm,BBOrd=BBOrd,indicI=indicI,Tmax=Tmax,orderTime=orderTime))
}



# PTCMestimMY
#############
# Fits a promotion time cure model with or without mismeasured covariates, using the backfitting approach of Ma and Yin (2008)
# Inputs:
  # - x: design matrix (of dimensions n x (R+P)) with no column of 1's
  # - y: matrix (of dimensions n x 2) containing all the ordered observed times in the 1st column, and the status of each observation in the 2nd column (0=censored, 1=failure, 2=cured)
  # - OmegaMatEstim: variance-covariance matrix (of dimensions (R+P+1)x(R+P+1)) of the measurement error 
  # - initMY: vector (of length R+P+1) containing initial values for the regression coefficients beta
  # - nBack: maximum number of iterations in the backfitting approach
  # - eps: convergence criterion
  # - multTimeMYestim: when multiplied by 300, gives the maximal number of seconds during which function nleqslv is allowed to run before being stopped.
  
PTCMestimMY <- function(x,y,OmegaMatEstim,initMY,nBack,eps,multTimeMYestim)
{
	enableJIT(3)

	Dat <- cbind(y,x)
	n <- dim(y)[1]
	P <- sum(diag(OmegaMatEstim)!=0)
	R <- dim(x)[2]-P # without error
	nbParam <- R+P
	variableNames <- colnames(x)

	varianceMY <- rep(NA,nbParam+1)

	prepData <- dataPreparation(Dat,P,R)
	OObsTime <- prepData$OObsTime
	V <- prepData$V
	Xc <- prepData$Xc
	mm <- prepData$mm
	BBOrd <- prepData$BBOrd
	indicI <- prepData$indicI
	Tmax <- prepData$Tmax

	MYres <- MaYinEstimation(mm,R,P,BBOrd,OObsTime,Xc,V,OmegaMatEstim=OmegaMatEstim,init=initMY,nBack=nBack,eps=eps,multTimeMYestim=multTimeMYestim,n=n)
	
	endK2 <- MYres$endK
	estimMY <- MYres$estimMYI
	mp <- MYres$mp
	flag <- MYres$flag

	# Variance estimation of the MY estimator	
	if(flag>0.5 & flag<1.5)
	{
		p <- rep(0,n)
		p[OObsTime[,2]==1] <- rep(mp,table(OObsTime[OObsTime[,2]==1,1]))
		fy <- matrix(nrow=n,ncol=mm,byrow=TRUE,data=(rep(BBOrd,n)<=rep(OObsTime[,1],rep(mm,n))))%*%mp
		fy[(OObsTime[,2]>1.5)] <- 1 #not necessary
	
		varianceMY <- computeVariance(x=cbind(1,V,Xc),bet=estimMY,sigma=OmegaMatEstim,d=OObsTime[,2],y=OObsTime[,1],fy=fy,n=n,nbParam=nbParam)
		colnames(varianceMY) <- rownames(varianceMY) <- c("Intercept",colnames(x))
	}
	return(list(coefficients=estimMY,estimCDF=fy,vcov=varianceMY,classObs=indicI,flag=flag,endK=endK2))
}



# End of: Functions called by PTCMestimBF

###############
# PTCMestimBF #
###############

PTCMestimBF <- function(x, ...) UseMethod("PTCMestimBF")


# PTCMestimBF.formula
#####################
# Fits a promotion time cure model with or without mismeasured covariates, using the backfitting approach of Ma and Yin (2008)
# Inputs:
  # - formula: 
  # - data: a dataframe containing, in the first 2 columns, the observed times and censoring indicators, and then the covariates
  # - ...: additional arguments to be passed to PTCMestimBF.default
  
PTCMestimBF.formula <- function(formula, data=list(), ...)
{
	mf <- model.frame(formula=formula, data=data)
	x <- model.matrix(attr(mf, "terms"), data=mf)
	x <- x[,-1]
	Y <- model.extract(mf, "response")
    if (!inherits(Y, "Surv")) 
        stop("Response must be a survival object")
	
	est <- PTCMestimBF.default(x, Y, ...)
	est$call <- match.call()
	est$formula <- formula
	est
}

# PTCMestimBF.default
#####################
# Fits a promotion time cure model with or without mismeasured covariates, using the backfitting approach of Ma and Yin (2008)
# Inputs:
  # - x: design matrix (of dimensions n x (R+P)) with no column of 1's
  # - y: matrix (of dimensions n x 2) containing all the ordered observed times in the 1st column, and the status of each observation in the 2nd column (0=censored, 1=failure, 2=cured)
  # - varCov: variance-covariance matrix (of dimensions (R+P+1)x(R+P+1)) of the measurement error 
  # - init: vector (of length R+P+1) containing initial values for the regression coefficients beta
  # - nBack: maximum number of iterations in the backfitting approach
  # - eps: convergence criterion
  # - multMaxTime: when multiplied by 300, gives the maximal number of seconds during which function nleqslv is allowed to run before being stopped.
  
PTCMestimBF.default <- function(x, y, varCov, init, nBack=10000, eps=1e-8, multMaxTime=2,...)
{
	x <- as.matrix(x)
	if (!inherits(y, "Surv")) 
        stop("Response must be a survival object")
	est <- PTCMestimMY(x, y, OmegaMatEstim=varCov, initMY=init, nBack=nBack, eps=eps, multTimeMYestim=multMaxTime)
	est$call <- match.call()
	class(est) <- "PTCMestimBF"
	est
}

# print.PTCMestimBF
###################
# Print an object of class "PTCMestimBF" to the screen
# Inputs: 
  # - x: an object of class PTCMestimBF, the result of function PTCMestimBF
  
print.PTCMestimBF <- function(x, ...)
{
	cat("Call:\n")
	print(x$call)
	cat("\nCoefficients:\n")
	print(x$coefficients)
}

# summary.PTCMestimBF
#####################
# Summarize a PTCMestimBF fit
# Inputs:
  # - object: an object of class PTCMestimBF, the result of function PTCMestimBF

summary.PTCMestimBF <- function(object, ...)
{
	se <- sqrt(diag(object$vcov))
	tval <- coef(object) / se
	TAB <- cbind(Estimate = coef(object),	StdErr = se,	z.value = tval,	p.value = 2*pnorm(-abs(tval)))
	res <- list(call=object$call,	coefficients=TAB)
	class(res) <- "summary.PTCMestimBF"
	res
}

# print.summary.PTCMestimBF
###########################
# Print the result of summary.PTCMestimBF
# Inputs:
  # - x: an object of class "summary.PTCMestimBF", the result of function summary.PTCMestimBF
  
print.summary.PTCMestimBF <- function(x, ...)
{
	cat("Call:\n")
	print(x$call)
	cat("\n")
	printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
}


######################################
# Functions called by PTCMestimSIMEX #
######################################

# SimexEstimation
#################
# Performs the first part (the Simulations part) of the SIMEX algorithm applied to the promotion time cure model
# Inputs:
  # - Nu: vector of levels of added noise
  # - B: number of replications for each level of added noise
  # - V: matrix containing the mismeasured covariates
  # - Xc: matrix containing the covariates without measurement error
  # - n: number of observations
  # - R: number of covariates without measurement error   
  # - P: number of mismeasured covariates
  # - OObsTime: matrix (of dimensions n x 2) containing all the ordered observed times in the 1st column, and the status of each observation in the 2nd column (0=censored, 1=failure, 2=cured)
  # - errorDistEstim: the distribution of the measurement error
  # - paramDistEstim: a scalar or a vector of length 2 containing the parameter(s) of the measurement error distribution, for nongaussian distributions
  # - OmegaMatEstim: variance-covariance matrix (of dimensions (R+P+1)x(R+P+1)) of the measurement error  
  # - mm: the number of distinct failure times
  # - BBOrd: vector (of length m) containing the ordered distinct failure times
  # - nBack: maximum number of iterations in the backfitting approach
  # - eps: convergence criterion
  # - init: vector (of length R+P+1) containing initial values for the regression coefficients beta
  # - multTimeMYestim: when multiplied by 300, gives the maximal number of seconds during which function nleqslv is allowed to run before being stopped. 
  
SimexEstimation <- function(Nu,B,V,Xc,n,R,P,OObsTime,errorDistEstim,paramDistEstim,OmegaMatEstim,mm,BBOrd,nBack,eps,init,multTimeMYestim)
{
	nbParam <- R+P
	indexNu <- 0
	estimNuIBF <- matrix(nrow=length(Nu),ncol=R+P+1)	
	estimBNuIBF <- matrix(nrow=B,ncol=(1+R+P)) 	
	varHatNuI <- matrix(nrow=length(Nu),ncol=R+P+1)
	varHatNuI[1,] <- 0
	varBarNuI <- matrix(nrow=length(Nu),ncol=R+P+1)
	varianceSimexBNu <- array(dim=c(B,R+P+1,R+P+1),data=NA)
	
	UU <- array(dim=c(B,n,P))
	
	if(length(paramDistEstim)>1)
	{
		paramDistEstim1 <- paramDistEstim[1]
		paramDistEstim2 <- paramDistEstim[2]
	}else
	{
		paramDistEstim1 <- paramDistEstim
	}

	for (nu in Nu) # for each level of added noise
	{
		indexNu <- indexNu+1
		print(paste("Level of added noise : ",nu))
		flush.console()

		if((nu-0)<0.000001) # When no added measurement error, no need to run B Simex-simulations: we obtain the naive estimator
		{
			VbNu1 <- V # nxp matrix, rows: indiv, cols: covariates
			WW <- cbind(1,array(dim=c(n,P),VbNu1),Xc)

			resBF <- MaYinEstimation(mm,R,P,BBOrd,OObsTime,Xc,VbNu1,OmegaMatEstim=matrix(nrow=R+P+1,ncol=R+P+1,data=0),init=init,nBack,eps,multTimeMYestim,n=n)
			estimNuIBF[indexNu,] <- resBF$estimMYI
			mp <- resBF$mp
		
			p <- rep(0,n)
			p[OObsTime[,2]==1] <- rep(mp,table(OObsTime[OObsTime[,2]==1,1]))
			fy <- matrix(nrow=n,ncol=mm,byrow=TRUE,data=(rep(BBOrd,n)<=rep(OObsTime[,1],rep(mm,n))))%*%mp
			fy[(OObsTime[,2]>1.5)] <- 1 #not necessary
			varBarNuI[indexNu,] <- diag(computeVariance(x=WW,bet=estimNuIBF[indexNu,],sigma=matrix(nrow=1+nbParam,ncol=1+nbParam,0),d=OObsTime[,2],y=OObsTime[,1],fy=fy,n=n,nbParam=nbParam))
		} # End of (if nu==0)
		
		if((nu-0)>=0.000001)
		{
			for (b in 1:B) # for each simulation (in one level of added noise)
			{
				# In the first OUTER loop (first non-zero level of added noise), add the noise U_{b,i} to the observed mismeasured covariates to obtain V_{b,i}(nu) for i=1,...,n (since U_{b,i} and V_{b,i} do not depend on nu, the level of added noise)
				if (nu==Nu[2]) 
				{
					UU[b,,] <- switch(errorDistEstim,normal=mvrnorm(n,mu=rep(0,P), Sigma=OmegaMatEstim[2:(1+P),2:(1+P)]),
							student=array(rt(n*P,paramDistEstim1)*paramDistEstim2,dim=c(n,P)),
							chiSquare=array((rchisq(n*P,paramDistEstim1)-paramDistEstim1)*paramDistEstim2,dim=c(n,P)),
							laplace=array(r(DExp(rate=1/paramDistEstim1))(n*P),dim=c(n,P))) 					
				}
					VbNub <- t(t(V)+sqrt(nu)%*%t(UU[b,,])) # The b-th matrix of Vb is an nxp matrix, rows: indiv, cols: covariates
	
				resBF <- MaYinEstimation(mm,R,P,BBOrd,OObsTime,Xc,V=array(dim=c(n,P),VbNub),OmegaMatEstim=matrix(nrow=R+P+1,ncol=R+P+1,data=0),init=estimNuIBF[1,],nBack,eps,multTimeMYestim,n=n)
				estimBNuIBF[b,] <- resBF$estimMYI
				mp <- resBF$mp
			
				p <- rep(0,n)
				p[OObsTime[,2]==1] <- rep(mp,table(OObsTime[OObsTime[,2]==1,1]))
				fy <- matrix(nrow=n,ncol=mm,byrow=TRUE,data=(rep(BBOrd,n)<=rep(OObsTime[,1],rep(mm,n))))%*%mp
				fy[(OObsTime[,2]>1.5)] <- 1 #not necessary

				varianceSimexBNu[b,,] <- computeVariance(x=cbind(1,V,Xc),bet=estimBNuIBF[b,],sigma=matrix(nrow=1+nbParam,ncol=1+nbParam,0),d=OObsTime[,2],y=OObsTime[,1],fy=fy,n=n,nbParam=nbParam)
			} # End of the B Simex-simulations
		
			# Compute the mean, over all simulations for this level of added noise, of each estimated parameter.
		
			estimNuIBF[indexNu,] <- apply(estimBNuIBF,2,mean) 
			varHatNuI[indexNu,] <- ((B-1)/B)*apply(estimBNuIBF,2,var)
			varBarNuI[indexNu,] <- diag(apply(varianceSimexBNu,c(2,3),mean))
		}
	}
	return(list(estimNuIBF=estimNuIBF, varHatNuI=varHatNuI, varBarNuI=varBarNuI))
}

# SimexExtrapolation
####################
# Performs the second part (the Extrapolation part) of the SIMEX algorithm applied to the promotion time cure model
# Inputs:
# - estimNui: a matrix containing the naive estimates for each parameter (the columns), for each value of added noise (the rows)
# - variableNames: vector of names for the covariates (including the intercept)
# - nbParam: total number of regression parameters
# - Nu: vector of levels of added noise
# - orderExtrap: order of the extrapolation function

SimexExtrapolation <- function(estimNui,variableNames,nbParam,Nu,orderExtrap)
{  
	estim <- rep(NA,1+nbParam)
	estimNuExMY <- matrix(nrow=length(Nu),ncol=orderExtrap+nbParam+1)
	estimNuExMY[,1:(nbParam+1)] <- estimNui
	for(j in 1:orderExtrap)
	{
		estimNuExMY[,j+nbParam+1] <- Nu^j
	}
	temp1 <- paste("nu",1:orderExtrap,sep="")
	temp2 <- paste(temp1,collapse="+")
	temp3 <- paste("~",temp2,collapse="")
	colnames(estimNuExMY) <- c(variableNames,temp1)
	temp4 <- (-1)^(1+1:(orderExtrap+1))

	for (p in 1:(1+nbParam)) # For each parameter to be estimated
	{
		fm2 <- as.formula(paste(colnames(estimNuExMY)[p],temp3))
		res <- lm(formula=fm2, data=data.frame(estimNuExMY))
		estim[p] <- temp4%*%res$coeff # Simex estimator 
	}
	return(list(estim=estim))
}


# PTCMestimSimex
################
# Fits a promotion time cure model with mismeasured covariates using the SIMEX method
# Inputs:
  # - x: design matrix (of dimensions n x (R+P)) with no column of 1's
  # - y: matrix (of dimensions n x 2) containing all the ordered observed times in the 1st column, and the status of each observation in the 2nd column (0=censored, 1=failure, 2=cured)
  # - errorDistEstim: the distribution of the measurement error
  # - paramDistEstim: a scalar or a vector of length 2 containing the parameter(s) of the measurement error distribution, for nongaussian distributions
  # - OmegaMatEstim: variance-covariance matrix (of dimensions (R+P+1)x(R+P+1)) of the measurement error 
  # - nBack: maximum number of iterations in the backfitting approach
  # - eps: convergence criterion
  # - Nu: vector of levels of added noise
  # - B: number of replications for each level of added noise
  # - initSimex: vector (of length R+P+1) containing initial values for the regression coefficients beta
  # - orderExtrap: a scalar or a numeric vector containing the degrees of the polynomials used in the extrapolation step
  # - multTimeMYestim: when multiplied by 300, gives the maximal number of seconds during which function nleqslv is allowed to run before being stopped.	
		
PTCMestimSimex <- function(x,y,errorDistEstim,paramDistEstim,OmegaMatEstim,nBack,eps,Nu,B,initSimex,orderExtrap,multTimeMYestim)
{
	enableJIT(3)

	Dat <- cbind(y,x)
	n <- dim(Dat)[1]
	P <- 1 # with error
	R <- dim(Dat)[2]-2-P # without error
	nbParam <- R+P
	variableNames <- c("Intercept",colnames(x))

	multTimeMYestim <- 2

	estimNuBF <- matrix(nrow=length(Nu),ncol=nbParam+1)
	estimSimexBF <- matrix(nrow=length(orderExtrap),ncol=1+nbParam,NA)
	colnames(estimNuBF) <- colnames(estimSimexBF) <- c("Intercept",colnames(x))

	varianceSimexBNu <- array(dim=c(length(Nu),B,nbParam+1,nbParam+1))
	varianceSimexBF <- matrix(nrow=length(orderExtrap),ncol=nbParam+1,NA)
	colnames(varianceSimexBF) <- c("Intercept",colnames(x))
	varHatNu <- matrix(nrow=length(Nu),ncol=nbParam+1)
	varHatNu[1,] <- 0
	varBarNu <- matrix(nrow=length(Nu),ncol=nbParam+1)

	prepData <- dataPreparation(Dat,P,R)
	OObsTime <- prepData$OObsTime
	V <- prepData$V
	Xc <- prepData$Xc
	mm <- prepData$mm
	BBOrd <- prepData$BBOrd
	indicI <- prepData$indicI
	Tmax <- prepData$Tmax
	
	resSimexMY <- SimexEstimation(Nu,B,V=V,Xc=Xc,n,R,P,OObsTime,errorDistEstim,paramDistEstim=paramDistEstim,OmegaMatEstim, mm,BBOrd,nBack,eps,init=initSimex,multTimeMYestim)

	estimNuBF <- resSimexMY$estimNuIBF
	varBarNu <- resSimexMY$varBarNuI
	varHatNu <- resSimexMY$varHatNuI	

	for(j in 1:length(orderExtrap))
	{
		quadExtrapMY <- SimexExtrapolation(estimNui=estimNuBF,variableNames=variableNames,nbParam=nbParam,Nu,orderExtrap[j])
		estimSimexBF[j,] <- quadExtrapMY$estim		

		varSimex <- SimexExtrapolation(estimNui=varBarNu-varHatNu,variableNames=variableNames,nbParam=nbParam,Nu,orderExtrap[j])
		varianceSimexBF[j,] <- varSimex$estim
	}
	return(list(coefficients=estimSimexBF,var=varianceSimexBF,classObs=indicI,estimNuBF=estimNuBF))
}

# End of: Functions called by PTCMestimSIMEX

##################
# PTCMestimSIMEX #
##################

PTCMestimSIMEX <- function(x, ...) UseMethod("PTCMestimSIMEX")

# PTCMestimSIMEX.formula
########################
# Fits a promotion time cure model with mismeasured covariates, using the SIMEX method
# Inputs:
  # - formula: 
  # - data: a dataframe containing, in the first 2 columns, the observed times and censoring indicators, and then the covariates
  # - ...: additional arguments to be passed to PTCMestimSIMEX.default

PTCMestimSIMEX.formula <- function(formula, data=list(), ...)
{
	mf <- model.frame(formula=formula, data=data)
	x <- model.matrix(attr(mf, "terms"), data=mf)
	x <- x[,-1]
	Y <- model.extract(mf, "response")
    if (!inherits(Y, "Surv")) 
        stop("Response must be a survival object")
	
	est <- PTCMestimSIMEX.default(x, Y, ...)
	est$call <- match.call()
	est$formula <- formula
	est
}

# PTCMestimSIMEX.default
########################
# Fits a promotion time cure model with or without mismeasured covariates, using the SIMEX method
# Inputs:
  # - x: design matrix (of dimensions n x (R+P)) with no column of 1's
  # - y: matrix (of dimensions n x 2) containing all the ordered observed times in the 1st column, and the status of each observation in the 2nd column (0=censored, 1=failure, 2=cured)
  # - errorDistEstim: the distribution of the measurement error
  # - paramDistEstim: a scalar or a vector of length 2 containing the parameter(s) of the measurement error distribution, for nongaussian distributions
  # - varCov: variance-covariance matrix (of dimensions (R+P+1)x(R+P+1)) of the measurement error 
  # - nBack: maximum number of iterations in the backfitting approach
  # - eps: convergence criterion
  # - Nu: vector of levels of added noise
  # - B: number of replications for each level of added noise
  # - init: vector (of length R+P+1) containing initial values for the regression coefficients beta
  # - orderExtrap: a scalar or a numeric vector containing the degrees of the polynomials used in the extrapolation step
  # - multMaxTime: when multiplied by 300, gives the maximal number of seconds during which function nleqslv is allowed to run before being stopped.
 
PTCMestimSIMEX.default <- function(x, y, errorDistEstim=c("normal","student","chiSquare","laplace"), paramDistEstim=NA, varCov=NA,  nBack=10000, eps=1e-8, Nu=c(0,.5,1,1.5,2), B=50, init, orderExtrap=2, multMaxTime=2,...)
{
	x <- as.matrix(x)
	if (!inherits(y, "Surv")) 
        stop("Response must be a survival object")	
	est <- PTCMestimSimex(x, y, errorDistEstim,paramDistEstim,varCov,nBack,eps,Nu,B,initSimex=init,orderExtrap,multTimeMYestim=multMaxTime)

	est$call <- match.call()
	est$orderExtrap <- orderExtrap
	class(est) <- "PTCMestimSIMEX"
	est
}

# print.PTCMestimSIMEX
######################
# Print an object of class "PTCMestimSIMEX" to the screen
# Inputs: 
  # - x: an object of class PTCMestimSIMEX, the result of function PTCMestimSIMEX
  
print.PTCMestimSIMEX <- function(x, ...)
{
	cat("Call:\n")
	print(x$call)
	for(i in 1:length(x$orderExtrap))
	{
		cat(paste("\nCoefficients with extrapolant of order ",x$orderExtrap[i],":\n"))
		print(x$coefficients[i,])
	}
}

# summary.PTCMestimSIMEX
########################
# Summarize a PTCMestimSIMEX fit
# Inputs:
  # - object: an object of class PTCMestimSIMEX, the result of function PTCMestimSIMEX
  
summary.PTCMestimSIMEX <- function(object, ...)
{
	res <- list(call=object$call,orderExtrap=object$orderExtrap)
	for(i in 1:length(object$orderExtrap))
	{
		se <- sqrt(object$var[i,])
		tval <- object$coefficients[i,] / se
		TAB <- cbind(Estimate = object$coefficients[i,],	StdErr = se,	z.value = tval,	p.value = 2*pnorm(-abs(tval)))
		res[[i+2]] <- list(coefficients=TAB)
	}
	class(res) <- "summary.PTCMestimSIMEX"
	res
	
}

# print.summary.PTCMestimSIMEX
##############################
# Print the result of summary.PTCMestimSIMEX
# Inputs:
  # - x: an object of class "summary.PTCMestimSIMEX", the result of function summary.PTCMestimSIMEX
  
print.summary.PTCMestimSIMEX <- function(x, ...)
{
	cat("Call:\n")
	print(x$call)
	for(i in 1:length(x$orderExtrap))
	{	
		cat(paste("\nResults with extrapolant of order ",x$orderExtrap[i],":\n"))
		printCoefmat(coef(x[[i+2]]), P.values=TRUE, has.Pvalue=TRUE)
	}
}
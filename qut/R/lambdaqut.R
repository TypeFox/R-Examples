lambdaqut <-
function(y,X,family=gaussian,alpha.level='default',M=1000,qut.standardize=TRUE,intercept=TRUE,no.penalty=NULL,offset=NULL,bootstrap=TRUE,beta0=NA){
#FUNCTION TO OBTAIN THE QUANTILE UNIVERSAL THRESHOLD

	#Hidden function for vectorized bootstrapping
	qutbootstrap <- function(curz){
		ind=sample(1:n,n,replace=TRUE)
		curX=X[ind,]
		curmuhat=curz[1:n]
		curz=curz[-(1:n)]
		if(!intercept&sum(O!=0)!=0&is.null(no.penalty)){
			curmuhat=muhatZ[ind,1]
		}
		else if(!is.null(no.penalty)|sum(O!=0)!=0){
			curbeta0=glm.fit(y=curz,x=as.matrix(A[ind,]),intercept=FALSE,family=family,offset=O[ind])$coef
			curmuhat=family$linkinv(as.matrix(A[ind,])%*%curbeta0+O[ind]) 
		}
		curbp=abs(t(curX)%*%(curz-curmuhat))
		
	}
	
	#Hidden function for vectorized computation of the intercept for each of the Monte Carlo simulations of the NULL model
	glm0=function(y,x,family,offset=offset){
		out=glm.fit(y=y,x=x,intercept=FALSE,family=family,offset=offset)
		return(out$coefficients)
	}
	
	family=family()
	n=nrow(X)
	p=ncol(X)
	X0=X
	if(is.null(offset)) O=rep(0,n)
	else O=offset

	#Check for warnings
	if(alpha.level!='default'){
		if(alpha.level>1|alpha.level<0){
			warning("alpha.level is not in [0,1] interval; set to default")
			alpha.level='default'
		}
	}
	if(M<=0){
		warning("M is <=0; set to 1000")
        M=1000
	}
	if (is.null(p) | (p <= 1)) stop("X should be a matrix with 2 or more columns")
	if (n!=length(y)) stop("Number of observations in y not equal to number of rows in X")
	if(length(O)!=n) stop("length of offset is different to the number of observations")
	if(family$family=='poisson'&(sum(y==0)==n)) stop("All your Poisson counts are zero")
	if(family$family=='binomial'&((sum(y==0)==n)|(sum(y==1)==n))) stop("All your Binomial measures are zero or one")
	
	#initialize A matrix
	if(!intercept&is.null(no.penalty)) muhat=family$linkinv(O)
	else{
		A=c()
		if(!is.null(no.penalty)){ 
			A=as.matrix(X[,no.penalty]) #if there are more unpenalized coefficients
			X=X[,-no.penalty]
		}
		if(intercept) A=cBind(rep(1,n),A)  #if there is an intercept (column of ones)

		#Estimate beta0 as glm(y~A)
		if(family$family=='gaussian') beta0=rep(0,ncol(A))
		else if(!is.numeric(beta0)) beta0=glm.fit(y=y,x=as.matrix(A),intercept=FALSE,family=family,offset=offset)$coef
		else if(length(beta0)!=(length(no.penalty)+intercept)) stop("length of beta0 is different from the number of unpenalized covariates or the intercept has not been included")
		muhat=family$linkinv(as.matrix(A)%*%beta0+O)
	}
		
	#Set default alpha.level
	if(alpha.level=='default')	alpha.level=1/(sqrt(pi*log(p)))

	#Check if A=1 and there is an explicit characterization of D
	if(!(intercept&is.null(no.penalty)&sum(O==0)==n)) warning("Explicit characterization of D is not defined")

	#Monte Carlo Simulation of the null model
	znotinD=0
	if(family$family=='gaussian') z=matrix(rnorm(n*M,mean=as.vector(muhat),sd=1),n,M)
	else if(family$family=='poisson'){
		z=matrix(rpois(n*M,lambda=as.vector(muhat)),n,M)
		znotinD=sum(apply(z,2,sum)==0)
		z=z[,apply(z,2,sum)!=0]
	}	
	else if(family$family=='binomial'){
		z=matrix(rbinom(n*M,size=1,prob=as.vector(muhat)),n,M)
		znotinD=sum((apply(z,2,sum)==0)|(apply(z,2,sum)==n))
		z=z[,(apply(z,2,sum)!=0)&(apply(z,2,sum)!=n)]
	}
	else stop("Not available family")

	if(znotinD>(M-2)) stop("Can't generate valid simulations under the null hypothesis, try increasing M")
	
	#obtain estimate of beta0 for each simulation
	if(!intercept&is.null(no.penalty)){
		muhatZ=family$linkinv(O)
		muhatZ=matrix(rep(muhatZ,ncol(z)),nrow=n)
	}
	else{
		beta0Z=apply(z,2,glm0,x=as.matrix(A),family=family,offset=offset)
		muhatZ=family$linkinv(as.matrix(A)%*%beta0Z+O)
	}

	divX=rep(1,ncol(X))
	
	if(!bootstrap){
		#No Bootstrap Matrix X
		
		#X'z
		bp=abs(t(X)%*%(z-muhatZ))
		scale.factor=rep(1,p)
		
		#quantile-based standardization
		if(qut.standardize){
			divX=apply(abs(bp),1,quantile,prob=(1-alpha.level))
			divX[which(divX==0)]=1
			
			#Alignment/scaling
			X=t(t(X)/divX)
			
			#X'z for the standardized X matrix
			bp=abs(t(X)%*%(z-muhatZ))
		}
	}
	else{
		#Bootstrap Matrix X

		#quantile-based standardization
		if(qut.standardize){
			bp=abs(t(X)%*%(z-muhatZ))
			scale.factor=rep(1,p)
			divX=apply(abs(bp),1,quantile,prob=(1-alpha.level))
			divX[which(divX==0)]=1
			
			#Alignment/scaling
			X=t(t(X)/divX)
		}

		#Bootstrapping
		bp=apply(rbind(muhatZ,z),2,qutbootstrap)
	}
	
	#Obtain distribution of Lambda=max(X'z)
	resultsMC=apply(bp,2,max,na.rm=TRUE)
	resultsMC=c(resultsMC,rep(Inf,znotinD))
	#obtain lambda.qut for the given quantile
	lambda=quantile(resultsMC, prob=(1-alpha.level))
	
	lambdamax=max(abs(t(X)%*%(y-muhat)))
	
	if(!is.null(no.penalty)){
		X0[,-no.penalty]=X
		scale.factor[-no.penalty]=divX
		X=X0
	}
	else scale.factor=divX

	#OUTPUT
	out=NULL
	out$scale.factor=scale.factor
	out$lambda.max=lambdamax
	out$lambda=lambda
	out$Xnew=X
	return(out)
	
}

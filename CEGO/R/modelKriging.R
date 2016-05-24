#   Copyright (c) 2014-2015 by Martin Zaefferer, Cologne University of Applied Sciences

###################################################################################
#' Kriging Model
#' 
#' Implementation of a distance-based Kriging model, e.g., for mixed or combinatorial input spaces.
#' It is based on employing suitable distance measures for the samples in input space.
#'
#' The basic Kriging implementation is based on the work of Forrester et al. (2008). 
#' For adaptation of Kriging to mixed or combinatorial spaces, as well as
#' choosing distance measures with Maximum Likelihood Estimation, see the other two references (Zaefferer et al., 2014).
#'
#' @param x list of samples in input space
#' @param y column vector of observations for each sample
#' @param distanceFunction a suitable distance function of type f(x1,x2), returning a scalar distance value, preferably between 0 and 1.
#'      Maximum distances larger 1 are no problem, but may yield scaling bias when different measures are compared.
#' 		Should be non-negative and symmetric.  It can also be a list of several distance functions. In this case, Maximum Likelihood Estimation (MLE) is used 
#'		to determine the most suited distance measure.
#' @param control (list), with the options for the model building procedure:\cr
#' \code{lower} lower boundary for theta, default is \code{1e-6}\cr
#' \code{upper} upper boundary for theta, default is \code{100}\cr
#' \code{corr} function to be used for correlation modelling, default is \code{fcorrGauss}\cr
#' \code{algTheta}  algorithm used to find theta (as well as p and lambda), default is \code{\link{optimInterface}}.\cr
#' \code{algThetaControl}  list of controls passed to \code{algTheta}.\cr
#' \code{optimizeP} boolean that specifies whether the exponents (\code{p}) should be optimized. Else they will be set to two. \cr
#' \code{useLambda} whether or not to use the regularization constant lambda (nugget effect). Default is \code{FALSE}.\cr
#' \code{lambdaLower} lower boundary for lambda, default is \code{-6}\cr 
#' \code{lambdaUpper} upper boundary for lambda, default is \code{0}\cr
#' \code{distances} a distance matrix. If available, this matrix is used for model building, instead of calculating the distance matrix using the parameters \code{distanceFunction}. Default is \code{NULL}.
#'
#' @return an object of class \code{modelKriging} containing the options and determined parameters for the model:\cr
#' \code{x} (see parameters)\cr
#' \code{y} (see parameters)\cr
#' \code{lower} (see parameters)\cr
#' \code{upper} (see parameters)\cr
#' \code{algTheta} (see parameters)\cr
#' \code{algThetaControl} (see parameters)\cr
#' \code{optimizeP} (see parameters)\cr
#' \code{theta} activity or width parameter theta, a parameter of the correlation function determined with MLE\cr
#' \code{log10Theta} log10 \code{theta} (i.e. \code{log10(theta)})\cr
#' \code{lambda} regularization constant (nugget) lambda \cr
#' \code{log10Lambda} log10 of regularization constant (nugget) lambda (i.e. \code{log10(lambda)})\cr
#' \code{p} exponent p, parameter of the correlation function determined with MLE (if \code{optimizeP} is \code{TRUE})\cr
#' \code{yMu} vector of observations y, minus MLE of mu\cr
#' \code{SSQ} Maximum Likelihood Estimate (MLE) of model parameter sigma^2\cr
#' \code{mu} MLE of model parameter mu\cr
#' \code{Psi} correlation matrix Psi\cr
#' \code{Psinv} inverse of Psi\cr
#' \code{nevals} number of Likelihood evaluations during MLE of theta/lambda/p\cr
#' \code{scaling} Default is FALSE. If TRUE: Distances values are divided by maximum distance to independent of the scale of the distance function.\cr
#' \code{distanceFunctionIndexMLE} If a list of several distance measures (\code{distanceFunction}) was given, this parameter contains the index value of the measure chosen with MLE.
#' 
#' @seealso \code{\link{predict.modelKriging}} 
#' 
#' @references Forrester, Alexander I.J.; Sobester, Andras; Keane, Andy J. (2008). Engineering Design via Surrogate Modelling - A Practical Guide. John Wiley & Sons.
#' @references Zaefferer, Martin; Stork, Joerg; Friese, Martina; Fischbach, Andreas; Naujoks, Boris; Bartz-Beielstein, Thomas. (2014). Efficient global optimization for combinatorial problems. In Proceedings of the 2014 conference on Genetic and evolutionary computation (GECCO '14). ACM, New York, NY, USA, 871-878. DOI=10.1145/2576768.2598282 http://doi.acm.org/10.1145/2576768.2598282 
#' @references Zaefferer, Martin; Stork, Joerg; Bartz-Beielstein, Thomas. (2014). Distance Measures for Permutations in Combinatorial Efficient Global Optimization. In Parallel Problem Solving from Nature - PPSN XIII (p. 373-383). Springer International Publishing.
#'
#' @examples
#' #Set random number generator seed
#' set.seed(1)
#' #Simple test landscape
#' fn <- landscapeGeneratorUNI(1:5,distancePermutationHamming)
#' #Generate data for training and test
#' x <- unique(replicate(40,sample(5),FALSE))
#' xtest <- x[-(1:15)]
#' x <- x[1:15]
#' #Determin true objective function values
#' y <- fn(x)
#' ytest <- fn(xtest)
#' #Build model
#' fit <- modelKriging(x,y,distancePermutationHamming,
#'     control=list(algThetaControl=list(method="L-BFGS-B"),useLambda=FALSE))
#' #Predicted obj. function values
#' ypred <- predict(fit,xtest)$y
#' #Uncertainty estimate
#' fit$predAll <- TRUE
#' spred <- predict(fit,xtest)$s
#' #Plot
#' plot(ytest,ypred,xlab="true value",ylab="predicted value",
#'     pch=20,xlim=c(0.3,1),ylim=c(min(ypred)-0.1,max(ypred)+0.1))
#' segments(ytest, ypred-spred,ytest, ypred+spred)
#' epsilon = 0.02
#' segments(ytest-epsilon,ypred-spred,ytest+epsilon,ypred-spred)
#' segments(ytest-epsilon,ypred+spred,ytest+epsilon,ypred+spred)
#' abline(0,1,lty=2)
#' #Use a different/custom optimizer (here: SANN) for maximum likelihood estimation: 
#' #(Note: Bound constraints are recommended, to avoid Inf values.
#' #This is really just a demonstration. SANN does not respect bound constraints.)
#' optimizer1 <- function(x,fun,lower=NULL,upper=NULL,control=NULL,...){
#'   res <- optim(x,fun,method="SANN",control=list(maxit=100),...)
#'   list(xbest=res$par,ybest=res$value,count=res$counts)
#' }
#' fit <- modelKriging(x,y,distancePermutationHamming,
#'                    control=list(algTheta=optimizer1,useLambda=FALSE))
#' #One-dimensional optimizer (Brent). Note, that Brent will not work when 
#' #several parameters have to be set, e.g., when using nugget effect (lambda).
#' #However, Brent may be quite efficient otherwise.
#' optimizer2 <- function(x,fun,lower,upper,control=NULL,...){
#'  res <- optim(x,fun,method="Brent",lower=lower,upper=upper,...)
#'  list(xbest=res$par,ybest=res$value,count=res$counts)
#' }
#' fit <- modelKriging(x,y,distancePermutationHamming,
#'                     control=list(algTheta=optimizer2,useLambda=FALSE))
#' @export
###################################################################################
modelKriging <- function(x, y, distanceFunction,control=list()){ #TODO: scaling for distance values not in [0;1]
	if(!is.matrix(y))
		y <- as.matrix(y)
	if(any(duplicated(x))){ #duplicates HAVE to be removed, but TODO: duplicates for noisy problems are okay. needs additional if(){} in relation to useLambda parameter?
		duplicates <- which(duplicated(x))
		x <- x[-duplicates]
		y <- as.matrix(y[-duplicates])
	}
		
	con<-list(lower=1e-6, upper=1e5, 
						corr=fcorrGauss, 
						algTheta= optimInterface, 
						algThetaControl= list(funEvals=200,reltol=1e-6,restarts=TRUE,method="L-BFGS-B"),#TODO: set by default to Brent, if one-dimensional?
						optimizeP= FALSE, 
						useLambda=FALSE, lambdaLower = -6, lambdaUpper = 0, 
						scaling=FALSE);
	con$algThetaControl[names(control$algThetaControl)] <- control$algThetaControl
	control$algThetaControl <- con$algThetaControl
	con[names(control)] <- control;
	control<-con;
	
	algThetaControl <- control$algThetaControl
	useLambda <- control$useLambda
	lambdaLower <- control$lambdaLower
	lambdaUpper <- control$lambdaUpper
	fcorr <- control$corr
  
	fit <- control
	fit$x <- x
	fit$y <- y

	lowerTheta <- log10(fit$lower)

	upperTheta <- log10(fit$upper)

	#Wrapper for optimizing theta  based on forrRegLikelihood:
	fitFun <- function (x, fX, fy,optimizeP,useLambda,fcorr){ 
		as.numeric(modelKrigingLikelihood(x,fX,fy,optimizeP,useLambda,fcorr)$NegLnLike)
		#print(c(as.numeric(result$NegLnLike),x))
		#return(as.numeric(result$NegLnLike)) #without return should be faster
	}

	n <- length(fit$x) #number of observations
	
	x0 <- fit$startTheta	

	#calculate distance matrix? or multiple matrices?
	if(length(distanceFunction)==1){
		if(is.null(control$distances))
			A <-distanceMatrix(x,distanceFunction) 
		else
			A <- control$distances
    maxD <- max(A) #maximum distance
    if(control$scaling){
      A <- A/maxD
    }    
	}else{
		if(is.null(control$distances)){
			A <- list()
      maxD <- list()
			for(i in 1:length(distanceFunction)){
        A[[i]] <-distanceMatrix(x,distanceFunction[[i]]) 
        maxD[[i]] <- max(A[[i]]) #maximum distance
        if(control$scaling){
          A[[i]] <- A[[i]]/maxD[[i]]
        } 
			}
		}else{
			A <- control$distances
      maxD <- list()
			for(i in 1:length(distanceFunction)){
        maxD[[i]] <- max(A[[i]]) #maximum distance
        if(control$scaling){
          A[[i]] <- A[[i]]/maxD[[i]]
        } 
			}    
		}
	}

	# start point for theta:	
	if(is.null(fit$startTheta))
		x1 <- n/100 #todo: works, but not sure whether good in general.
	else
		x1 <- fit$startTheta

	#todo: it has to be noted, that in original kriging the design space is scaled, thus also the distance matrix.
	# this yields always distances in the same range, and theta values in the same range (if they have a similar importance).
	# Here, no scaling is performed. The drawback is, that this may yield the requirement of completely different 
	# theta value ranges for different problems. I.e.: The lower and upper boundary for theta is hard to specify.
	# Thus, it may be advisable to scale the yielded distances to a range of 0 to 1.
	#A <- A/max(A) #scaling with maximum difference? need to rescale in predictor
	
	if(fit$optimizeP){ # optimize p
		lowerTheta <- c(lowerTheta, 0.01)
		upperTheta <- c(upperTheta, 2)		
		x3 <- 1 #start values for p
		x0 <- c(x1,x3)
	}else{ # p  is fixed to 1 
		x0 <- c(x1)
	}
	if(useLambda){
		# start value for lambda:
		x2 <- lambdaLower + (lambdaUpper - lambdaLower)*runif(1)
		x0 <- c(x0,x2)
		#append regression constant lambda (nugget)
		lowerTheta <- c(lowerTheta,lambdaLower)
		upperTheta <- c(upperTheta, lambdaUpper)
	}	

	#force x0 into bounds
	x0 <- pmin(x0,upperTheta)
	x0 <- pmax(x0,lowerTheta)
	algThetaControl$funEvals <- algThetaControl$funEvals*length(x0)
	
	if(length(distanceFunction)==1){
		res <- control$algTheta(x=x0,fun=fitFun,lower=lowerTheta,upper=upperTheta,
						control=algThetaControl,fX=A,fy=fit$y,optimizeP=fit$optimizeP,useLambda=useLambda,fcorr=fcorr)	
		fit$distanceFunction <- distanceFunction
	}else{
		res <- list()
		minlik=Inf
		minlikindex=1
		for(i in 1:length(distanceFunction)){
			res[[i]] <- control$algTheta(x=x0,fun=fitFun,lower=lowerTheta,upper=upperTheta,
						control=algThetaControl,fX=A[[i]],fy=fit$y,optimizeP=fit$optimizeP,useLambda=useLambda,fcorr=fcorr)	
			if(res[[i]]$ybest < minlik){
				minlik <- res[[i]]$ybest
				minlikindex <- i
			}
		}	
		res <- res[[minlikindex]]
		fit$distanceFunction <- distanceFunction[[minlikindex]]
		A<-A[[minlikindex]]
    maxD <- maxD[[minlikindex]]
		fit$distanceFunctionIndexMLE <- minlikindex
	}	
	if(is.null(res$xbest))res$xbest <- x0;
	Params <- res$xbest
	nevals <- as.numeric(res$count[[1]])
	# extract model parameters:
	fit$theta <- 10^Params[1]
	if(fit$optimizeP){	
		fit$p <- Params[2]
	}	
	fit$log10Theta <- Params[1]
	if(useLambda){
		fit$log10Lambda <- Params[length(Params)];
		fit$lambda <- 10^Params[length(Params)]
	}else{
		fit$log10Lambda <- -Inf;
		fit$lambda <- 0
	}
	res <- modelKrigingLikelihood(c(fit$log10Theta,fit$p, fit$log10Lambda),A,fit$y,fit$optimizeP,useLambda,fcorr)
	
	fit$yMu <- res$yMu
	fit$SSQ <- as.numeric(res$SSQ)
	fit$mu <- res$mu
	fit$Psi <- res$Psi
	fit$Psinv <- res$Psinv
	fit$nevals <- nevals
	fit$like <- res$NegLnLike
  fit$maximumDistance <- maxD
  fit$predAll <- FALSE
	class(fit)<- "modelKriging"
	if(is.na(fit$Psinv[1])){ #model building failed. no invertible correlation matrix was found. return NA fit
		stop("Building the Kriging model failed, no invertible correlation matrix was found. This may be due to the specific data-set or distance function used.")
	}else{
		return(fit)
	}	
}

###################################################################################
#' Kriging Model
#'
#' DEPRECATED version of the Kriging model, please use \code{\link{modelKriging}}
#' 
#' @param x list of samples in input space
#' @param y column vector of observations for each sample
#' @param distanceFunction a suitable distance function of type f(x1,x2), returning a scalar distance value
#' @param control options for the model building procedure
#'
#' @keywords internal
#' @export
###################################################################################
combinatorialKriging <- function(x, y, distanceFunction, control = list()){
	.Deprecated("modelKriging")
	modelKriging(x,y,distanceFunction,control)
}

###################################################################################
#' Calculate negative log-likelihood
#' 
#' Used to determine theta/lambda/p values for the Kriging model in \code{\link{modelKriging}}
#' with Maximum Likelihood Estimation (MLE).
#'
#' @param xt vector, containing log10(theta), p and lambda
#' @param dX matrix of distances/dissimilarites between training samples
#' @param y vector of observations at sample locations
#' @param optimizeP whether to optimize p or not (FALSE at default)
#' @param useLambda whether to use nugget effect, i.e., lambda (FALSE at default)
#' @param corr whether to use nugget effect, i.e., lambda (fcorrGauss at default)
#'
#' @return list with elements\cr
#' \code{NegLnLike}  concentrated log-likelihood *-1 for minimising \cr
#' \code{Psi} correlation matrix\cr
#' \code{Psinv} inverse of correlation matrix (to save computation time in forrRegPredictor)\cr
#' \code{mu} MLE of model parameter mu \cr
#' \code{yMu} vector of observations y minus mu\cr
#' \code{SSQ} MLE of model parameter sigma^2
#'
#' @seealso \code{\link{modelKriging}}
#' @keywords internal
###################################################################################
modelKrigingLikelihood <- function(xt,dX,y,optimizeP=FALSE,useLambda=FALSE,corr=fcorrGauss){
	#browser()
	if(optimizeP){
		dX <- abs(dX)^(xt[2])
	}
	theta <- 10^xt[1];
	n <- dim(y)[1] #number of observations	
	Psi <- corr(theta,dX)
	if(useLambda){
		lambda <- 10^xt[length(xt)];
		Psi <- Psi + diag(lambda,n) 
	}		
	
	## cholesky decomposition
	cholPsi <- try(chol(Psi), TRUE) 
	
	## use pivoting if standard fails
	if(class(cholPsi) == "try-error"){
		cholPsi <- try(chol(Psi,pivot=TRUE), TRUE) 
	}	
	
	## give penalty if both fail
	if(class(cholPsi) == "try-error"){
		warning("Correlation matrix is not positive semi-definite (in modelKrigingLikelihood). Returning penalty.")
		return(list(NegLnLike=1e4,Psi=NA,Psinv=NA,mu=NA,SSQ=NA))
	}	
		
	#calculate natural log of the determinant of Psi (numerically more reliable and also faster than using det or determinant)
	LnDetPsi <- 2*sum(log(abs(diag(cholPsi))))
	
	#inverse with cholesky decomposed Psi
	Psinv <- try(chol2inv(cholPsi),TRUE)
	if(class(Psinv) == "try-error"){
		warning("Correlation matrix is not positive semi-definite (in modelKrigingLikelihood). Returning penalty.")
		return(list(NegLnLike=1e4,Psi=NA,Psinv=NA,mu=NA,SSQ=NA))
	}	
	
	
	psisum <- sum(Psinv) #this sum of all matrix elements may sometimes become zero, which may be caused by inaccuracies. then, the following may help
	if(psisum==0){
		psisum <- as.numeric(rep(1,n) %*% Psinv %*% rep(1,n))
		if(psisum==0){ #if it is still zero, return penalty
			warning("Sum of elements in inverse correlation matrix is zero (in modelKrigingLikelihood). Returning penalty.")
			return(list(NegLnLike=1e4,Psi=NA,Psinv=NA,mu=NA,SSQ=NA))
		}
	}		
	mu <- sum(Psinv%*%y)/psisum
	yMu <- y-mu 
	SigmaSqr <- (t(yMu)%*%Psinv%*%yMu)/n
	if(SigmaSqr < 0){ #TODO: or better return absolute?
		warning("Maximum Likelihood Estimate of model parameter sigma^2 is negative (in modelKrigingLikelihood). Returning penalty. ")
		return(list(NegLnLike=1e4-SigmaSqr,Psi=NA,Psinv=NA,mu=NA,SSQ=NA)) 
	}
	NegLnLike <- n*log(SigmaSqr) + LnDetPsi
	if(is.infinite(NegLnLike))#this may happen if all y are 0
		return(list(NegLnLike=1e4,Psi=NA,Psinv=NA,mu=NA,SSQ=NA)) 
	list(NegLnLike=NegLnLike,Psi=Psi,Psinv=Psinv,mu=mu,yMu=yMu,SSQ=SigmaSqr)
}


###################################################################################
#' Predict: Combinatorial Kriging
#' 
#' Predict with a model fit resulting from \code{\link{modelKriging}}.
#'
#' @param object fit of the Kriging model (settings and parameters), of class \code{modelKriging}.
#' @param x list of samples to be predicted
#' @param ... further arguments, not used
#'
#' @return Returned value depends on the setting of \code{object$predAll}\cr
#' TRUE: list with function value (mean) \code{$y} and uncertainty estimate \code{$s} (standard deviation)\cr
#' FALSE:\code{$y}only
#'
#' @seealso \code{\link{modelKriging}}
#' @export
###################################################################################
predict.modelKriging <- function(object,x,...){ #TODO: re-interpolation
	if(!is.list(x))x<-list(x)
	xo <- object$x
	theta <- object$theta 
	Psinv <- object$Psinv 
	n <- length(xo)
	#one <- rep(1,n)
	mu <- object$mu
	yMu <- object$yMu	
	SigmaSqr <- object$SSQ
	psi <- matrix(1,length(x),n)
	fundist <-object$distanceFunction
	if(object$optimizeP){
		p <- object$p
		for (i in 1:n)
			psi[,i] <- distanceVector(xo[[i]],x,fundist)^p
	}else{	
		for (i in 1:n)
			psi[,i] <- distanceVector(xo[[i]],x,fundist)#unlist(lapply(x,fndist))#^2
		###or in one line, but this is slower.
		#psi <- t(outer(xo,x,function(a,b)mapply(fundist,a,b))^2)
	}	
  if(object$scaling){
    psi <- psi/object$maximumDistance
  }
	psi <- object$corr(theta,psi)
	##
	f <- as.numeric(psi%*%(Psinv%*%yMu))+mu 
	## return value:
	res <- list(y=f)
	##########################################################################
	if (object$predAll){
		lambda <- object$lambda
		SSqr <- SigmaSqr*(1+lambda-diag(psi%*%(Psinv%*%t(psi)))) 
		s <- sqrt(abs(SSqr))
		res$s <- as.numeric(s)
	}
	res
}

###################################################################################
#' Gaussian Kernel for Kriging
#'
#' @param theta kernel parameter
#' @param D distance matrix
#'
#' @return matrix (Psi)
#'
#' @seealso \code{\link{modelKriging}}
#' 
#' @export
#' @keywords internal
###################################################################################
fcorrGauss <- function(theta,D){
	exp(-theta * D)
}

###################################################################################
#' Cubic Kernel for Kriging
#'
#' @param theta kernel parameter
#' @param D distance matrix
#'
#' @return matrix (Psi)
#'
#' @seealso \code{\link{modelKriging}}
#' 
#' @export
#' @keywords internal
###################################################################################
fcorrCubic <- function(theta,D){
	Psi = pmin(D * theta,1)
	1 - Psi^2 * (3 - 2*Psi)
}

###################################################################################
#' Linear Kernel for Kriging
#'
#' @param theta kernel parameter
#' @param D distance matrix
#'
#' @return matrix (Psi)
#'
#' @seealso \code{\link{modelKriging}}
#' 
#' @export
#' @keywords internal
###################################################################################
fcorrLinear <- function(theta,D){
	pmax(1- D * theta,0)
}

###################################################################################
#' Spherical Kernel for Kriging
#'
#' @param theta kernel parameter
#' @param D distance matrix
#'
#' @return matrix (Psi)
#'
#' @seealso \code{\link{modelKriging}}
#' 
#' @export
#' @keywords internal
###################################################################################
fcorrSphere <- function(theta,D){
	Psi = pmin(D * theta,1)
	 1 - Psi * (1.5 - 0.5*Psi^2)
}	

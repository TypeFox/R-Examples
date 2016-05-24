###################################################################################
#' Meta Model Interface: Forrester's Kriging
#' 
#' Interface to the Kriging model based on Matlab code by Forrester et al. 2008.
#' Can be used both for single and multi objective SPOT.
#'
#' @param rawB unmerged data
#' @param mergedB merged data
#' @param design new design points which should be predicted
#' @param spotConfig global list of all options, needed to provide data for calling functions. The model specific settings in this list are\cr
#'	\code{spotConfig$seq.forr.loval} lower boundary for theta, default is \code{1e-3}\cr
#'	\code{spotConfig$seq.forr.upval} upper boundary for theta, default is \code{100}\cr
#'	\code{spotConfig$seq.forr.opt.p} boolean that specifies whether the exponents (\code{p}) should be optimized. Else they will be set to two. Default value is \code{FALSE}. Default is highly recommended as the implementation of this feature is not yet well tested and might be faulty.\cr
#'	\code{spotConfig$seq.forr.algtheta} algorithm used to find theta, default is \code{"optim-L-BFGS-B"}. Else, any from the list of possible \code{method} values in \code{\link{spotOptimizationInterface}} can be chosen.\cr
#'	\code{spotConfig$seq.forr.budgetalgtheta} budget for the above mentioned algorithm, default is \code{100}. The value will be multiplied with the length of the model parameter vector to be optimized.
#'	\code{spotConfig$seq.forr.reinterpolate} boolean that specifies whether re-interpolation should be used during the prediction process. Default value is \code{FALSE}. Setting this to \code{TRUE} is recommended, when an error estimate of nearly zero is desired at sample locations, regardless of chosen regularization constant (nugget). Please note that prediction with interpolation will take longer than without.\cr
#'	\code{spotConfig$seq.forr.savetheta} boolean that specifies whether the exponents (\code{p}) should be optimized. Else they will be set to two. Default value is \code{FALSE}. Default is recommended since this feature not yet well tested, and might lead to a preference of local optima.\cr
#' @param fit if an existing model fit is supplied, the model will not be build based on 
#'				data, but only evaluated with the model fit (on the design data). To build the model, 
#'				this parameter has to be NULL. If it is not NULL the parameters mergedB and rawB will not be 
#'				used at all in the function.
#'
#' @return returns the list \code{spotConfig} with two new entries:\cr
#' 	spotConfig$seq.modelFit fit of the model used with the predictor functions \cr
#'	spotConfig$seq.largeDesignY the y values of the design, evaluated with the fit
#' @export
#' @seealso \code{\link{forrBuilder}} \code{\link{forrRegPredictor}} \code{\link{forrReintPredictor}} 
#' @references Forrester, Alexander I.J.; Sobester, Andras; Keane, Andy J. (2008). Engineering Design via Surrogate Modelling - A Practical Guide. John Wiley & Sons.
###################################################################################
spotPredictForrester <- function(rawB,mergedB,design,spotConfig,fit=NULL){	
	design <- spotInitializePredictor(design,"data.frame",spotConfig$alg.roi,NULL,"spotPredictForrester",spotConfig$io.verbosity)	#MZ: removed MASS. not necessary?
	if(is.null(spotConfig$seq.forr.reinterpolate))spotConfig$seq.forr.reinterpolate=FALSE
	########################################################
	# BUILD
	########################################################	
	if(is.null(fit)){
		if(is.null(spotConfig$seq.forr.control))spotConfig$seq.forr.control=list()
		
		xNames <- row.names(spotConfig$alg.roi)
		yNames <- spotConfig$alg.resultColumn
		x <- as.matrix(rawB[xNames])
		if(length(yNames)==1){
			y <- as.matrix(rawB[[yNames]])
			listXY <- spotRepairMissingValues(x,y,spotConfig)
			y = as.matrix(listXY$y)
			x = listXY$x
			spotConfig$seq.na.fit <- listXY$nafit
			fit <- forrBuilder(x, y, spotConfig$alg.roi$lower, spotConfig$alg.roi$upper, spotConfig$seq.forr.control)
		}
		else{#Distinction for multi criteria spot 
			yMat <- rawB[yNames]
			fit=list()
			spotConfig$seq.na.fit <- list()
			for (i in 1:length(yNames)){
				y<-as.matrix(yMat[,i])
				listXY <- spotRepairMissingValues(x,y,spotConfig)
				y = as.matrix(listXY$y)
				x = listXY$x
				spotConfig$seq.na.fit[[i]] <- listXY$nafit				
				fit[[i]]<-forrBuilder(x, y,spotConfig$alg.roi$lower, spotConfig$alg.roi$upper, spotConfig$seq.forr.control)
				}			
		}
	}
	########################################################
	# PREDICT
	########################################################
	if(!is.null(design)){ 		
		pred.all=spotConfig$seq.model.variance
		nmodel <- length(spotConfig$alg.resultColumn)
		if(nmodel>1){ #do multi criteria prediction
			resy=matrix(0,nrow(design),nmodel)
			resvar=matrix(NA,nrow(design),nmodel)
			y=list()
			for (i in 1:length(fit)){ #predict			
				if(spotConfig$seq.forr.reinterpolate){
					res= forrReintPredictor(design,fit[[i]],pred.all)
				}else{
					res= predict(fit[[i]],design,pred.all)
				}
				if(!is.null(spotConfig$seq.na.penalty)) # assign a penalty
					resy[,i]= spotPenalizeMissingValues(design,res$f,spotConfig$seq.na.fit[[i]], spotConfig$seq.na.penalty[i])
				else #no penalty
					resy[,i]= res$f
				if(pred.all)resvar[,i]= res$s
				y[[i]]= fit[[i]]$y
			}
			if(is.function(spotConfig$seq.infill)){# do EI 
				resy= spotConfig$seq.infill(resy,resvar,y,spotConfig$mco.refPoint)
			}
		}else{ #do single criteria prediction
			if(spotConfig$seq.forr.reinterpolate){
				res= forrReintPredictor(design,fit,pred.all)
			}else{
				res <- predict(fit,design,pred.all)
			}
			if(!is.null(spotConfig$seq.na.penalty)) # assign a penalty
				resy= spotPenalizeMissingValues(design,res$f,spotConfig$seq.na.fit, spotConfig$seq.na.penalty)
			else # no penalty
				resy= res$f			
			resvar=matrix(NA,nrow(design),1)
			if(pred.all)resvar=res$s
			if(is.function(spotConfig$seq.infill)){ # do EI			
				resy= spotConfig$seq.infill(resy,resvar,min(fit$y))
			}
		}			
	}else{
		resy <- NULL
		resvar <- NULL
	}	
	########################################################
	# OUTPUT
	########################################################	
	spotWriteLines(spotConfig$io.verbosity,3,"spotPredictForrester finished")
	spotConfig$seq.largeDesignY=as.data.frame(resy)	
	spotConfig$seq.largeDesignVar=as.data.frame(resvar)	
	spotConfig$seq.modelFit<-fit;
	spotConfig
}

###################################################################################
#' Print Function Forrester Kriging
#'
#' Print information about a Forrester Kriging fit, as produced by \code{\link{forrBuilder}}.
#'
#' @rdname print
#' @method print forr
# @S3method print forr
#' @param x	fit returned by \code{\link{forrBuilder}}.
#' @param ... additional parameters	
#' @export
#' @keywords internal
###################################################################################
print.forr <- function(x,...){
	cat("------------------------\n")
	cat("Forrester Kriging model.\n")
	cat("------------------------\n")
	cat("Estimated activity parameters (theta) sorted \n")
	cat("from most to least important variable \n")
	cat(paste("x",order(x$dmodeltheta,decreasing=TRUE),sep="",collaps=" "))
	cat("\n")	
	cat(sort(x$dmodeltheta,decreasing=TRUE))	
	cat("\n \n")
	cat("exponent(s) p:\n")
	if(x$opt.p)
		cat(x$P)
	else
		cat(2)
	cat("\n \n")
	cat("Estimated regularization constant (or nugget) lambda:\n")	
	cat(x$dmodellambda)
	cat("\n \n")
	cat("Number of Likelihood evaluations during MLE:\n")	
	cat(x$nevals)
	cat("\n")	
	cat("------------------------\n")
}

###################################################################################
#' Build Forrester Kriging
#'
#' This function builds a Kriging model based on code by Forrester et al..
#' By default exponents (p) are fixed at a value of two, and a nugget (or regularization constant) is used.
#' 
#' @param X design matrix (sample locations)
#' @param y vector of observations at X
#' @param lb lower boundary of the design space. Will be extracted from the matrix \code{X} if not given.
#' @param ub upper boundary of the design space. Will be extracted from the matrix \code{X} if not given.
#' @param control (list), with the options for the model building procedure:\cr
#' \code{loval} lower boundary for theta, default is \code{1e-6}\cr
#' \code{upval} upper boundary for theta, default is \code{100}\cr
#' \code{corr} function to be used for correlation modeling, default is \code{fcorrGauss}\cr
#' \code{algtheta}  algorithm used to find theta, default is \code{"NLOPT_LN_NELDERMEAD"}. Else, any from the list of possible \code{method} values in \code{spotOptimizationInterface} from the SPOT package can be chosen.\cr
#' \code{budgetalgtheta} budget for the above mentioned algorithm, default is \code{100}. The value will be multiplied with the length of the model parameter vector to be optimized.
#' \code{opt.p} boolean that specifies whether the exponents (\code{p}) should be optimized. Else they will be set to two. \cr
#' \code{uselambda} whether or not to use the regularization constant lambda (nugget effect). Default is \code{FALSE}.
#' \code{lambda.loval} lower boundary for lambda, default is \code{-6}\cr 
#' \code{lambda.upval} upper boundary for lambda, default is \code{0}\cr
#' \code{starttheta} optional start value for theta.
#' \code{reinterpolate} whether (TRUE) or not (FALSE, default) reinterpolation should be performed
#'
#' @return a fit (list), with the options and found parameters for the model which has to be passed to the predictor function:\cr
#' \code{X} sample locations (scaled to values between 0 and 1)\cr
#' \code{y} observations at sample locations (see parameters)\cr
#' \code{loval} lower boundary for theta (see parameters)\cr
#' \code{upval} upper boundary for theta (see parameters)\cr
#' \code{algtheta} algorithm to find theta (see parameters)\cr
#' \code{budgetalgtheta} budget for the above mentioned algorithm (see parameters)\cr
#' \code{opt.p} boolean that specifies whether the exponents (\code{p}) were optimized (see parameters)\cr
#' \code{normalizeymin} minimum in normalized space\cr
#' \code{normalizeymax} maximum in normalized space\cr
#' \code{normalizexmin} minimum in input space\cr
#' \code{normalizexmax} maximum in input space\cr
#' \code{dmodeltheta} vector of activity parameters\cr
#' \code{Theta} log_10 vector of activity parameters (i.e. \code{log10(dmodeltheta)})\cr
#' \code{dmodellambda} regularization constant (nugget) \cr
#' \code{Lambda} log_10 of regularization constant (nugget) (i.e. \code{log10(dmodellambda)})\cr
#' \code{yonemu} \code{Ay-ones*mu} \cr
#' \code{ssq} sigma square\cr
#' \code{mu} mean mu\cr
#' \code{Psi} matrix large Psi\cr
#' \code{Psinv} inverse of Psi\cr
#' \code{nevals} number of Likelihood evaluations during MLE
#'
#' @export
#' @seealso \code{\link{spotPredictForrester}} \code{\link{predict.forr}} \code{\link{forrReintPredictor}}
#' @references Forrester, Alexander I.J.; Sobester, Andras; Keane, Andy J. (2008). Engineering Design via Surrogate Modelling - A Practical Guide. John Wiley & Sons.
#'
#' @examples
#' ## Create design points
#' x = cbind(runif(20)*15-5,runif(20)*15)
#' ## Compute observations at design points (for Branin function)
#' y = as.matrix(apply(x,1,spotBraninFunction))
#' ## Create model with default settings
#' fit = forrBuilder(x,y)
#' ## Print model parameters
#' print(fit)
#'
###################################################################################
forrBuilder <- function(X, y, lb=NULL, ub=NULL, control=list()){#,debug=NA){
	con<-list(loval=1e-3, upval=1e2, algtheta="optim-L-BFGS-B", budgetalgtheta=100, opt.p= FALSE, uselambda=TRUE, lambda.loval = -6, lambda.upval = 0, starttheta=NULL, reinterpolate=FALSE);
	con[(namc <- names(control))] <- control;
	control<-con;
	
	fit = control
	k = ncol(X)
	fit$X = X
	fit$y = y
	# normalize input data
	ymin = 0
	ymax = 1
	fit$normalizeymin=ymin
	fit$normalizeymax=ymax
	res = spotNormalizeMatrix(fit$X, ymin, ymax, lb, ub)
	fit$X =res$y
	fit$normalizexmin=res$xmin
	fit$normalizexmax=res$xmax
	LowerTheta = rep(1,k)*log10(fit$loval)

	UpperTheta = rep(1,k)*log10(fit$upval)

	#Wrapper for optimizing theta  based on forrRegLikelihood:
	fitFun <- function (x, fX, fy,opt.p,uselambda){ #todo vectorize, at least for cma_es with active vectorize?
		as.numeric(forrRegLikelihood(x,fX,fy,opt.p,uselambda)$NegLnLike)
	}
	n=nrow(fit$X) #number of observations
	
	
	if(is.null(fit$starttheta))
		x1 =  rep(n/(100*k),k) # start point for theta
	else
		x1 = fit$starttheta

	# start value for lambda:
	x2 = fit$lambda.loval + (fit$lambda.upval - fit$lambda.loval)*runif(1)   #todo also assign start values by user?
	
	A=matrix(0,k,n*n)

	for(i in 1:k){
		A[i,]=as.numeric(as.matrix(dist(fit$X[,i]))) #MZ: speedup fix, using dist function: 100%
	}
	
	if(fit$opt.p){ # optimize p
		LowerTheta = c(LowerTheta, rep(1,k)*0.01)
		UpperTheta = c(UpperTheta, rep(1,k)*2)		
		x3 = rep(1,k)* 1.9 #start values for p
		x0 = c(x1,x3)
	}else{ # p  is fixed to 2 and the array A is completely precalculated
		A=A^2
		x0 = c(x1)
	}			
	if(fit$uselambda){
		# start value for lambda:
		x2 = fit$lambda.loval + (fit$lambda.upval - fit$lambda.loval)*runif(1)
		x0 = c(x0,x2)
		#append regression constant lambda (nugget)
		LowerTheta = c(LowerTheta,fit$lambda.loval)
		UpperTheta = c(UpperTheta,fit$lambda.upval)
	}	
	#force x0 into bounds
	x0= pmin(x0,UpperTheta)
	x0= pmax(x0,LowerTheta)
	opts=list(fevals=fit$budgetalgtheta*length(x0), reltol=1e-6, restarts=TRUE)	
	#browser()
	res <- spotOptimizationInterface(par=x0,fn=fitFun,gr=NULL,lower=LowerTheta,upper=UpperTheta,method=fit$algtheta,
						control=opts,fX=A,fy=fit$y,opt.p=fit$opt.p,uselambda=fit$uselambda)	
	if(is.null(res$par))res$par=x0;
	Params = res$par
	nevals = as.numeric(res$counts[[1]])
	
	fit$dmodeltheta=10^Params[1:k]
	if(fit$opt.p){	
		fit$P=Params[(k+1):(2*k)]
	}
	if(fit$uselambda){
		fit$Lambda = Params[length(Params)];
		fit$dmodellambda=10^Params[length(Params)]
	}else{
		fit$Lambda = -Inf;
		fit$dmodellambda=0
	}
	# extract model parameters:
	fit$Theta = Params[1:k]
	res=forrRegLikelihood(c(fit$Theta,fit$P, fit$Lambda),A,fit$y,fit$opt.p,fit$uselambda);
	
	fit$yonemu=res$yonemu
	fit$ssq=as.numeric(res$ssq)
	fit$mu=res$mu
	fit$Psi=res$Psi
	fit$Psinv=res$Psinv
	fit$nevals=nevals
	fit$like=res$NegLnLike
	class(fit)<- "forr"
	fit
}

###################################################################################
#' Normalize design
#'
#' Normalize design with given maximum and minimum in input space
#' 
#' @param x design matrix in input space (n rows for each point, k columns for each parameter)
#' @param ymin minimum vector of normalized space
#' @param ymax maximum vector of normalized space
#' @param xmin minimum vector of input space
#' @param xmax maximum vector of input space
#'
#' @return normalized design matrix
#' @seealso \code{\link{spotPredictForrester}} \code{\link{forrBuilder}}
#'  \code{\link{spotNormalizeMatrix}} 
#' \code{\link{spotReverseNormalizeMatrix}}
#' @export
#' @keywords internal
###################################################################################
spotNormalizeMatrix2 <- function (x,ymin,ymax,xmin,xmax){ 
	rangex = xmax-xmin
	rangey = ymax-ymin
	s = dim(x)[1]
	rangey * (x-matrix(rep(xmin,s),nrow=s,byrow=TRUE))/matrix(rep(rangex,s),nrow=s,byrow=TRUE) + ymin
}

###################################################################################
#' Normalize design
#' 
#' Normalize design by using minimum and maximum of the design values for input space
#'
#' @param x design matrix in input space
#' @param ymin minimum vector of normalized space
#' @param ymax maximum vector of normalized space
#'
#' @return normalized design matrix
#' @seealso \code{\link{spotPredictForrester}} \code{\link{forrBuilder}}
#'  \code{\link{spotNormalizeMatrix2}}
#' \code{\link{spotReverseNormalizeMatrix}} 
#' @export
#' @keywords internal
###################################################################################
spotNormalizeMatrix <- function(x,ymin, ymax, xmin=NULL, xmax=NULL){
	# Return the maximum from each row:
	if(is.null(xmax))
		xmax <- apply(x,2,max)
	# Return the minimum from each row:
	if(is.null(xmin))
		xmin <- apply(x,2,min)
	s <- dim(x)[1]
	rangex <- xmax-xmin
	rangey <- ymax-ymin
	xmin[rangex==0] <- xmin[rangex==0]-0.5
	xmax[rangex==0] <- xmax[rangex==0]+0.5
	rangex[rangex==0] <- 1
	y <- rangey * (x-matrix(rep(xmin,s),nrow=s,byrow=TRUE))/matrix(rep(rangex,s),nrow=s,byrow=TRUE) + ymin
	list(y=y,xmin=xmin,xmax=xmax)
}

###################################################################################
#' Reverse Normalize
#' 
#' Reverse normalization of a design matrix
#'
#' @param y design matrix in norm space
#' @param ymin minimum vector of normalized space
#' @param ymax maximum vector of normalized space
#' @param xmin minimum vector of output space
#' @param xmax maximum vector of output space
#'
#' @return denormalized design matrix
#' @seealso \code{\link{spotPredictForrester}} \code{\link{forrBuilder}}
#'  \code{\link{spotNormalizeMatrix2}} \code{\link{spotNormalizeMatrix}} 
#' @export
#' @keywords internal
###################################################################################
spotReverseNormalizeMatrix <- function(y,xmin,xmax,ymin,ymax){
	rangex <- xmax-xmin
	rangey <- ymax-ymin
	rangex * (y-ymin)*(1/rangey) + xmin
}

###################################################################################
#' Calculate negative log-likelihood
#' 
#' Used to determine theta/lambda values for the Kriging model in \code{\link{forrBuilder}}.
#'
#' @param x vector, containing log10(theta) and lambda
#' @param AX 3 dimensional array, constructed by forrBuilder from the sample locations
#' @param Ay vector of observations at sample locations
#'
#' @return list with elements\cr
#' \code{NegLnLike}  concentrated log-likelihood *-1 for minimising \cr
#' \code{Psi} correlation matrix\cr
#' \code{Psinv} inverse of correlation matrix (to save computation time in predict.forr)\cr
#FIXMZ: \code{U} U matrix of the LU decomposition of correlation matrix \cr
#' \code{mu} \cr
#' \code{ssq}
#' @seealso \code{\link{spotPredictForrester}} \code{\link{forrBuilder}}
#' @export
#' @keywords internal
###################################################################################
forrRegLikelihood <- function(x,AX,Ay,opt.p=FALSE,uselambda=TRUE){
	nx<-nrow(AX)
	theta=10^x[1:nx]
	if(opt.p){	
		#	if(any(is.na(abs(AX)^(10^x[(nx+1):(2*nx)]))))stop("NA values in theta")
		AX=abs(AX)^(x[(nx+1):(2*nx)])
	}
	lambda=0
	if(uselambda)
		lambda=10^x[length(x)]
	if( any(theta==0) ||  any(is.infinite(c(theta,lambda)))){ #for instance caused by bound violation
		return(list(NegLnLike=1e4,Psi=NA,Psinv=NA,mu=NA,ssq=NA))
	}
	n=dim(Ay)[1]
	
	Psi=exp(-matrix(colSums(theta*AX),n,n))
	if(uselambda)
		Psi=Psi+diag(lambda,n)

	## cholesky decomposition
	cholPsi <- try(chol(Psi), TRUE) 
	
	## use pivoting if standard fails
	if(class(cholPsi) == "try-error"){
		cholPsi <- try(chol(Psi,pivot=TRUE), TRUE) 
	}	
	
	## give penalty if both fail
	if(class(cholPsi) == "try-error"){
		warning("Correlation matrix is not positive semi-definite (in combinatorialKrigingLikelihood). Returning penalty.")
		return(list(NegLnLike=1e4,Psi=NA,Psinv=NA,mu=NA,SSQ=NA))
	}	
		
	#calculate natural log of the determinant of Psi (numerically more reliable and also faster than using det or determinant)
	LnDetPsi=2*sum(log(abs(diag(cholPsi))))
	
	#inverse with cholesky decomposed Psi
	Psinv= try(chol2inv(cholPsi),TRUE)
	if(class(Psinv) == "try-error"){
		warning("Correlation matrix is not positive semi-definite (in combinatorialKrigingLikelihood). Returning penalty.")
		return(list(NegLnLike=1e4,Psi=NA,Psinv=NA,mu=NA,SSQ=NA))
	}		
	
	mu=sum(Psinv%*%Ay)/sum(Psinv)# note: matrix%*%onevector is signif. faster than rowSums(matrix)
	yonemu=Ay-mu 
	SigmaSqr=(t(yonemu)%*%Psinv%*%yonemu)/n
	NegLnLike=n*log(SigmaSqr) + LnDetPsi
	list(NegLnLike=NegLnLike,Psi=Psi,Psinv=Psinv,mu=mu,yonemu=yonemu,ssq=SigmaSqr)
}

###################################################################################
#' Predict Forrester Model
#' 
#' Predict samples on a Forrester Kriging model.
#'
#' @param object Kriging model (settings and parameters) of class forr
#' @param x design matrix to be predicted
#' @param predictAll if TRUE return all (RMSE and prediction, in a dataframe), else return only prediction
#' @param ... not used
#'
#' @return Returned value is dependent on the setting of \code{pred.all}\cr
#' TRUE: data.frame with columns f (function values) and s (RMSE)\cr
#' FALSE: vector of function values only
#'
#' @examples
#' ## Create design points
#' x = cbind(runif(20)*15-5,runif(20)*15)
#' ## Compute observations at design points (for Branin function)
#' y = as.matrix(apply(x,1,spotBraninFunction))
#' ## Create model
#' fit = forrBuilder(x,y)
#' ## Create candidate design points
#' xx = cbind(runif(20)*15-5,runif(20)*15)
#' ## Predict candidates
#' y1 = predict(fit,xx)$f
#' ## Plot model (in comments, due to time consumption)
#' #fn <- function(x){predict(fit,as.matrix(x))$f}
#' #spotSurf3d(fn,c(-5,0),c(10,15))
#' ## Plot real function
#' #spotSurf3d(function(x){apply(x,1,spotBraninFunction)},c(-5,0),c(10,15))
#'
#' @seealso \code{\link{forrBuilder}} \code{\link{forrReintPredictor}}
#' @export
###################################################################################
predict.forr <- function(object,x,predictAll=FALSE,...){
	if(object$reinterpolate){
		return(forrReintPredictor(x,object,predictAll))
	}
	#normalize input x
	x <- spotNormalizeMatrix2(as.matrix(x),0,1,object$normalizexmin,object$normalizexmax)
	AX=object$X
	#Ay=object$y
	theta=object$dmodeltheta
	#theta=rep(sum(theta),length(theta))
	Psinv=object$Psinv #fixed: does not need to be computed, is already done in likelihood function
	n=dim(AX)[1]
	#one=rep(1,n)
	mu=object$mu
	yonemu=object$yonemu	
	SigmaSqr=object$ssq
	psi=matrix(1,nrow(x),n);
	if(object$opt.p){
		p=object$P
		for (i in 1:n) #todo this is for each center, calculate distance to new samples. BUT this should be case sensitive: only forEachCenter if more samples than centers,else other way round?
			psi[,i]=colSums(theta*(abs(AX[i,]-t(x))^p))
	}else{
		for (i in 1:n)
			psi[,i]=colSums(theta*((AX[i,]-t(x))^2))
	}	
	psi = exp(-psi)
	f=as.numeric(psi%*%(Psinv%*%yonemu))+mu #vectorised
	##########################################################################
	#if (object$Option!="Pred"){
	if (predictAll){
		lambda=object$dmodellambda;
		SSqr= SigmaSqr*(1+lambda-diag(psi%*%(Psinv%*%t(psi)))) #vectorised
		#TODO "diag(psi%*%...)" is excessive, since diag wastes alot of values computed by %*%
		s=sqrt(abs(SSqr));
	}
	if(!predictAll){list(f=f)}else{data.frame(f=f,s=as.numeric(s))}
}

# for backwards compatibility:
###################################################################################
#' Predict Forrester Model
#' 
#' This function is for backwards compatibiltiy only. See new predictor: \code{\link{predict.forr}} 
#'
#' @param x design matrix to be predicted
#' @param fit Kriging model (settings and parameters) of class forr
#' @param pred.all if TRUE return all (RMSE and prediction, in a dataframe), else return only prediction
#'
#' @return Returned value is dependent on the setting of \code{pred.all}\cr
#' TRUE: data.frame with columns f (function values) and s (RMSE)\cr
#' FALSE: vector of function values only
#'
#' @seealso \code{\link{predict.forr}} 
#' @export
###################################################################################
forrRegPredictor<- function(x,fit,pred.all=FALSE){
	predict.forr(fit,x,pred.all)
} 

###################################################################################
#' Predict Forrester Model (Re-interpolating)
#' 
#' Kriging predictor with re-interpolation to avoid stalling the optimization process which employs this model as a surrogate.
#' This is supposed to be used with deterministic experiments, which do need a non-interpolating model that avoids predicting non-zero error at sample locations.
#' This can be useful when the model is deterministic (i.e. repeated evaluations of one parameter vector do not yield different values) but does have a "noisy" structure (e.g. due to computational inaccuracies, systematical error).
#'
#' Please note that this re-interpolation implementation will not necessarily yield values of exactly zero at the sample locations used for model building. Slight deviations can occur.
#'
#' @param x design matrix to be predicted
#' @param ModelInfo fit of the Kriging model (settings and parameters)
#' @param pred.all if TRUE return all (RMSE and prediction, in a dataframe), else return only prediction
#'
#' @return Returned value is dependent on the setting of \code{pred.all}\cr
#' TRUE: data.frame with columns f (function values) and s (RMSE)\cr
#' FALSE: vector of function values only
#'
#' @examples
#' ## Create design points
#' x = cbind(runif(20)*15-5,runif(20)*15)
#' ## Compute observations at design points (for Branin function)
#' y = as.matrix(apply(x,1,spotBraninFunction))
#' ## Create model
#' fit = forrBuilder(x,y)
#' ## first estimate error with regressive predictor
#' sreg = predict(fit,x,TRUE)$s
#' ## now estimate error with re-interpolating predictor
#' sreint = forrReintPredictor(x,fit,TRUE)$s
#' print(sreg)
#' print(sreint)
#' ## sreint should be close to zero, significantly smaller than sreg
#'
#' @seealso \code{\link{forrBuilder}} \code{\link{forrCoBuilder}} \code{\link{predict.forr}}
#' @export
###################################################################################
forrReintPredictor <- function(x,ModelInfo,pred.all=FALSE){
	#normalize input x
	x <- spotNormalizeMatrix2(as.matrix(x),0,1,ModelInfo$normalizexmin,ModelInfo$normalizexmax)
	AX=ModelInfo$X
	#Ay=ModelInfo$y
	theta=ModelInfo$dmodeltheta
	Psi=ModelInfo$Psi
	Psinv=ModelInfo$Psinv 
	lambda=ModelInfo$dmodellambda;
	n=dim(AX)[1]
	#one=rep(1,n)
	mu=ModelInfo$mu	
	yonemu=ModelInfo$yonemu	
	#
	PsiB=Psi-diag(lambda,n)+diag(.Machine$double.eps,n)
	SigmaSqr=(t(yonemu)%*%Psinv%*%PsiB%*%Psinv%*%yonemu)/n;
	#	
	psi=matrix(1,nrow(x),n);
	if(ModelInfo$opt.p){
		p=ModelInfo$P
		for (i in 1:n)
			psi[,i]=exp(-colSums(theta*(abs(AX[i,]-t(x))^p)))
	}else{
		for (i in 1:n)
			psi[,i]=colSums(theta*((AX[i,]-t(x))^2))
	}	
	psi = exp(-psi)
	f=as.numeric(psi%*%(Psinv%*%yonemu))+mu #vectorised
	##########################################################################
	#if (ModelInfo$Option!="Pred"){
	if(pred.all){
		#
		Psinv= try(solve.default(PsiB), TRUE) ##Important notes: chol2inv(chol(Psi)) may be less likely to give an answer BUT may also be faster and more accurate
		if(class(Psinv) == "try-error"){
			Psinv=ginv(Psi)
		}	
		#
		SSqr= SigmaSqr*(1-diag(psi%*%(Psinv%*%t(psi)))) #vectorised
		s=sqrt(abs(SSqr));
	}
	if(!pred.all){list(f=f)}else{data.frame(f=f,s=as.numeric(s))}
}


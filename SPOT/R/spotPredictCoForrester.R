#TODO: external functions directly as surrogate models -> write an interface


###################################################################################
#' Meta Model Interface: Forrester's Co-Kriging
#' 
#' Interface to the Co-Kriging model based on Matlab code by Forrester et al. 2008.
#' Make sure to define a cheap, correlated function to be used here in \code{spotConfig$seq.forr.co.fn}.
#' This should be in form of y=f(x). This function is not allowed to yield the same values as the actual (expensive) target function.
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
#' @seealso \code{\link{forrCoBuilder}} \code{\link{forrCoRegPredictor}} 
#' @references FORRESTER, A.I.J, SOBESTER A. & KEAN, A.J. (2007), Multi-Fidelity optimization via surrogate modelling. \emph{Proc. R. Soc. A} 463, 3251-3269. \cr
#' LE GRATIET, L. & GARNIER, J. (2012), Recursive co-kriging model for Design of Computer Experiments with multiple levels of fidelity, \emph{arXiv:1210.0686} \cr
#' Forrester, Alexander I.J.; Sobester, Andras; Keane, Andy J. (2008). Engineering Design via Surrogate Modelling - A Practical Guide. John Wiley & Sons.
###################################################################################
spotPredictCoForrester <- function(rawB,mergedB,design,spotConfig,fit=NULL) {
	design <- spotInitializePredictor(design,"data.frame",spotConfig$alg.roi,NULL,"spotPredictCoForrester",spotConfig$io.verbosity)	
	########################################################
	# BUILD
	########################################################	
	if(is.null(fit)){
		xNames <- row.names(spotConfig$alg.roi)
		yNames <- spotConfig$alg.resultColumn
		if(is.null(spotConfig$seq.forr.control))spotConfig$seq.forr.control=list()
		
			if(length(yNames)==1){	#TODO: fine function is single objective
			#  sample locations:
			xe <- as.matrix(spotConfig$alg.currentResult[xNames])	 #rawB not used here, because it is sorted by response value	
			#  observations:
			ye <- as.matrix(spotConfig$alg.currentResult[yNames])
			# cheap sample locations
			xc <- as.matrix(spotConfig$alg.currentCoResult[xNames])
			indf <- spotConfig$alg.currentCoResult$indf
			xc <- rbind(xc[!indf,,drop=FALSE],xc[indf,,drop=FALSE])
			# cheap observations: 
			yc <- as.matrix(spotConfig$alg.currentCoResult[yNames])
			yc <- rbind(yc[!indf,,drop=FALSE],yc[indf,,drop=FALSE])
			#repair yc and xc before building the coarse-only kriging model
			config <- spotConfig
			config$seq.modelFit <- spotConfig$seq.modelFit$c #here, only the cheap model fit can be used.
			config$seq.predictionModel.func <- "spotPredictForrester" #with the normal kriging model		
			listXY <- spotRepairMissingValues(xc,yc,config)
			yc1 <- as.matrix(listXY$y)
			xc1 <- listXY$x
			#
			if(is.null(spotConfig$seq.modelFit)|(is.null(spotConfig$seq.forr.contfitC))){		
				fitC <- forrBuilder(xc1, yc1, spotConfig$alg.roi$lower, spotConfig$alg.roi$upper, spotConfig$seq.forr.control)
			}else{ #if allready fit, and desired: use old theta value, and small budget for update
				controls <- spotConfig$seq.forr.control
				controls$starttheta <- spotConfig$seq.modelFit$c$Theta
				controls$budgetalgtheta <-  controls$budgetalgtheta / 5
				fitC <- forrBuilder(xc1, yc1, spotConfig$alg.roi$lower, spotConfig$alg.roi$upper, controls)		
			}			
			#
			listCE <- spotRepairMissingValuesCoKriging(xe, ye, xc, yc, spotConfig, fitC) #xc/yc have to be repaired now in a different way, because they can now be replaced even in the initial step by fitC
			yc <- as.matrix(listCE$yc) #todo these yc/xc values are not even used in the forrCoBuilder... if they are, there may be problems, because of differences with the yc/xc in fitC?
			xc <- listCE$xc			
			ye <- as.matrix(listCE$ye)
			xe <- listCE$xe		
			#note: it is expected that the cheap evaluations of the expensive sample locations are at the end of the xc/yc matrices
			#
			fit <- forrCoBuilder(xe, ye, xc, yc, fitC,  spotConfig$alg.roi$lower, spotConfig$alg.roi$upper, spotConfig$seq.forr.control);		
			#
		}else{NULL} #TODO else does not exist yet, MCO is a problem here. It has to be defined somehow, which objective values have a correlated variable from another target function (or whether there is correlation between themselves???)
	}
	########################################################
	# PREDICT
	########################################################
	if(!is.null(design)){ 		
		pred.all<-spotConfig$seq.model.variance
		nmodel <- length(spotConfig$alg.resultColumn)
		if(nmodel>1){ #do multi criteria prediction
			NULL #TODO
		}else{ #do single criteria prediction
			#design <- spotNormalizeMatrix2(t(as.matrix(design)),0,1,fit$normalizexmin,fit$normalizexmax);
			#if(spotConfig$seq.forr.reinterpolate){
			#	res<- forrCoReintPredictor(design,fit,pred.all)
			#}else{
				res<- forrCoRegPredictor(design,fit,pred.all)
			#}
			resy<-res$f
			resvar<-matrix(NA,nrow(design),1)
			if(pred.all)resvar<-res$s
			if(is.function(spotConfig$seq.infill)){ # do EI			
				resy<- spotConfig$seq.infill(resy,resvar,min(fit$y))
			}
		}
	}else{
		resy <- NULL
		resvar <- NULL
	}	
	########################################################
	# OUTPUT
	########################################################		
	spotConfig$seq.largeDesignY<-as.data.frame(resy)
	spotConfig$seq.largeDesignVar<-as.data.frame(resvar)	
	spotConfig$seq.modelFit<-fit;
	spotConfig
}

###################################################################################
#' Call Correlated Function
#'
#' This function is used to evaluate a cheap function (cheaper as the real target).
#' The purpose is to enhance a Kriging model with Co-Kriging.
#' 
#' @param spotConfig set of parameters, containing \code{spotConfig$seq.forr.co.fn}, which is the correlated (cheap) function.
#' @param result data frame of observations and related information
#' @param design design matrix to be evaluated on the correlated function.
#' @param dfunc co function to be evaluated
#'
#' @return normalized design matrix
#' @seealso \code{\link{spotPredictCoForrester}} \code{\link{forrBuilder}} \code{\link{forrCoBuilder}}
#' @keywords internal
###################################################################################
spotCallCoFunction <- function(spotConfig,result,design,dfunc){ #TODO remove, replace ... ?
	spotConfig$spot.fileMode <- FALSE
	spotConfig$alg.tar.func <- dfunc
	spotConfig$alg.func <- "spotOptimInterface"
	#####spotConfig$alg.resultColumn  ?
	spotConfig$io.apdFileName <- "noapdfile.apd"
	spotConfig$io.resFileName <- "noresfile.apd"
	spotConfig$alg.currentResult <- result
	spotConfig$alg.currentDesign <- design
	spotConfig$io.verbosity <- 0
	spotConfig$alg.seed <- NA #todo: should be setable?
	#####calculate predictions of the correlated function
	spotConfig<-do.call(spotConfig$alg.func, args=list(spotConfig))
	spotConfig$alg.currentResult
}

###################################################################################
#' helper function
#' @param x object
#' @param i index
#' @keywords internal
###################################################################################
forrIF <- function(x,i){
	ifelse(length(x)>1,x[[i]],x)
}

###################################################################################
#' Build Forrester Co-Kriging
#'
#' This function builds a Co-Kriging model based on code by Forrester et al..
#' Please note that the expensive sample locations should be contained in the cheap sample locations.
#' Furthermore, it has to be made sure that the correlated functions do not yield identical values.
#' That is, \code{ye} and \code{yc} should have common sample locations, but different values.
#' The sample locations only evaluated on the cheap function can be chosen arbitrarily.
#' 
#' @param Xe design matrix (expensive sample locations)
#' @param ye 1-row matrix of expensive observations at Xe
#' @param Xc design matrix (cheap sample locations). The bottom of this matrix should contain expensive samples.
#' @param yc 1-row matrix of cheap observations at Xc. 
#' @param fitC object of class \code{forr}, containing a Kriging model build through the cheap observations
#' @param lb lower boundary of the design space. Will be extracted from the matrix \code{Xe} if not given.
#' @param ub upper boundary of the design space. Will be extracted from the matrix \code{Xe} if not given.
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
#' \code{rho.loval} lower boundary for rho, default is \code{-5}\cr 
#' \code{rho.upval} upper boundary for rho, default is \code{5}\cr
#' \code{starttheta} optional start value for theta.
#' \code{reinterpolate} whether (TRUE) or not (FALSE, default) reinterpolation should be performed
#'
#' @return a fit (list) of class \code{coforr}. This contains Co-Kriging specific parameters, as well as two fits of class \code{forr} which represent the cheap and expensive models.
#'
#' @export
#' @seealso \code{\link{spotPredictCoForrester}} \code{\link{forrCoRegPredictor}} \code{\link{forrBuilder}}
#' @references Forrester, Alexander I.J.; Sobester, Andras; Keane, Andy J. (2008). Engineering Design via Surrogate Modelling - A Practical Guide. John Wiley & Sons.
#'
#' @examples
#'	## This is the one-variable example described by Forrester et al.
#'	## The "expensive" function to be modeled is
#'	ovar <- function(x){(x*6-2)^2*sin((x*6-2)*2)}
#'	## The "cheap" function to be modeled is
#'	covar <- function(x){ A=0.5;B=10;C=-5;D=0;
#'		A*(((x+D)*6-2)^2)*sin(((x+D)*6-2)*2)+((x+D)-0.5)*B+C
#'	}
#'	## construct cheap and expensive sample locations
#'	xe <- rbind(0,0.4,0.6,1)
#'	xc <- rbind(0.1,0.2,0.3,0.5,0.7,0.8,0.9,0,0.4,0.6,1)
#'	## get observations of samples
#'	ye <- rbind(ovar(xe))
#'	yc <- rbind(covar(xc))
#' 	## build the Co-Kriging model, with cheap and expensive observations
#'	set.seed(1)
#'  fitC <- forrBuilder(xc, yc)
#'	fit <- forrCoBuilder(xe, ye, xc, yc, fitC)
#' 	## build the ordinary Kriging model with expensive observations only
#'	fit1 <- forrBuilder(xe, ye)	 
#'  ## Predict and plot over whole design space
#'	x<-seq(from=0,to=1,by=0.01)
#'	yco <- forrCoRegPredictor(as.matrix(x),fit,FALSE)$f
#'	ypc <- predict(fitC,as.matrix(x),FALSE)$f
#' 	fit2 <- fit
#'  class(fit2) <- "forr"
#'	ype <- predict(fit2,as.matrix(x),FALSE)$f
#'	yy <- predict(fit1,as.matrix(x),FALSE)$f
#'	plot(x,ovar(x),type="l",ylim=c(-15,20),lwd=3)
#'	points(xe,ye,pch=19,cex=1.5)
#'	points(xc,yc,cex=1.5)
#'	lines(x,covar(x),lwd=3)	
#' 	lines(x,ype,col="blue",lwd=3) #difference model 
#'	lines(x,ypc,col="red",lty=4,lwd=3) 	#cheap model
#'	lines(x,yy,col="blue",lty=3,lwd=3)#uncorrected model
#'	lines(x,yco,col="darkgreen",lty=5,lwd=3) #comodel
#'	legend("top",lwd=c(3,3,1,1,3,3,3,3),
#'		col=c("black","black","black","black","blue","red","blue","darkgreen"),
#'		legend=	c("Expensive Function", "Cheap Function", 
#'				"Expensive Observations", "Cheap Observations",
#'				"Uncorrected Model", "Cheap Model","Difference Model",
#'				"Co-Kriging Model"),
#'		lty=c(1,1,0,0,3,4,1,5),pch=c(NA,NA,19,1,NA,NA,NA,NA))
#'	sum((yco-ovar(x))^2)/length(x) #mse
#'
#' @references FORRESTER, A.I.J, SOBESTER A. & KEAN, A.J. (2007), Multi-Fidelity optimization via surrogate modelling. \emph{Proc. R. Soc. A} 463, 3251-3269. \cr
#' Forrester, Alexander I.J.; Sobester, Andras; Keane, Andy J. (2008). Engineering Design via Surrogate Modelling - A Practical Guide. John Wiley & Sons.
###################################################################################
forrCoBuilder <- function(Xe,ye, Xc, yc, fitC, lb=NULL, ub=NULL, control=list()){
	
	con<-list(loval=1e-3, upval=1e2, algtheta="optim-L-BFGS-B", budgetalgtheta=100, opt.p= FALSE, uselambda=TRUE, lambda.loval = -6, lambda.upval = 0, starttheta=NULL, rho.loval=-5, rho.upval=5, reinterpolate=FALSE);
	con[(namc <- names(control))] <- control;
	control<-con;
	
	#########################################
	#first build a standard model for the cheap observations
	#########################################
	fitE= control
	k = ncol(Xe);
	fitE$X = Xe; 
	fitE$y = ye;
	# normalize input data
	ymin = 0;
	ymax = 1;
	fitE$normalizeymin=ymin;
	fitE$normalizeymax=ymax;
	res = spotNormalizeMatrix(fitE$X, ymin, ymax, lb, ub);
	fitE$X =res$y;
	fitE$normalizexmin=res$xmin;
	fitE$normalizexmax=res$xmax;

	LowerTheta = rep(1,k)*log10(fitE$loval);
	UpperTheta = rep(1,k)*log10(fitE$upval);

	#Wrapper for optimizing theta  based on forrRegLikelihood:
	fitFun <- function (x, fX, fy, ffy, opt.p,uselambda){ #todo vectorize, at least for cma_es with active vectorize?
		as.numeric(forrCoLikelihood(x,fX,fy,ffy,opt.p,uselambda)$NegLnLike)
	}
	n=nrow(fitE$X) #number of observations
	
	if(is.null(fitE$starttheta))
		x1 =  rep(n/(100*k),k) # start point for theta
	else
		x1 = fitE$starttheta

	#start value for rho 
	x4 = fitE$rho.loval + (fitE$rho.upval - fitE$rho.loval)*runif(1)
	
	
	
	####instead of the old solution above, concatenate matrices. is faster in likelihood function.
	A=matrix(0,k,n*n)
	for(i in 1:k){
		A[i,]=as.numeric(as.matrix(dist(fitE$X[,i]))) #MZ: speedup fix, using dist function: 100%
	}	
	if(fitE$opt.p){ # optimize p
		LowerTheta = c(LowerTheta, rep(1,k)*0.01)
		UpperTheta = c(UpperTheta, rep(1,k)*2)		
		x3 = rep(1,k)* 1.9 #start values for p
		x0 = c(x1,x3,x4)
	}else{ # p  is fixed to 2 and the array A is completely precalculated
		A=A^2
		x0 = c(x1,x4)
	}		
	#append rho bounds
	LowerTheta = c(LowerTheta, fitE$rho.loval);	
	UpperTheta = c(UpperTheta, fitE$rho.upval);		
	if(fitE$uselambda){
		# start value for lambda:
		x2 = fitE$lambda.loval + (fitE$lambda.upval - fitE$lambda.loval)*runif(1)
		x0 = c(x0,x2)
		#append regression constant lambda (nugget)
		LowerTheta = c(LowerTheta, fitE$lambda.loval);	
		UpperTheta = c(UpperTheta, fitE$lambda.upval);		
	}		

	lc=dim(fitC$y)[1] 

	opts=list(fevals=fitE$budgetalgtheta*length(x0), reltol=1e-6, restarts=TRUE)	
	res <- spotOptimizationInterface(par=x0,fn=fitFun,gr=NULL,lower=LowerTheta,upper=UpperTheta,method=fitE$algtheta,
						control=opts,fX=A,fy=fitE$y, ffy=fitC$y[(lc-n+1):lc], opt.p=fitE$opt.p,uselambda=fitE$uselambda) #todo what if fitC$y and yc are not the same?
	if(is.null(res$par))res$par=x0;
	Params = res$par
	nevals = as.numeric(res$counts[[1]])
	
	# extract model parameters:	
	fitE$dmodeltheta=10^Params[1:k];
	if(fitE$opt.p){	
		fitE$P=Params[(k+1):(2*k)];		
	}
	fitE$Theta = Params[1:k];
	if(fitE$uselambda){
		fitE$Lambda = Params[length(Params)]
		fitE$dmodellambda=10^Params[length(Params)]
		fitE$rho = Params[length(Params)-1]
		res=forrCoLikelihood(c(fitE$Theta,fitE$P,fitE$rho,fitE$Lambda),A,fitE$y,fitC$y[(lc-n+1):lc],fitE$opt.p,fitE$uselambda);
	}else{
		fitE$Lambda = -Inf
		fitE$dmodellambda=0
		fitE$rho = Params[length(Params)]
		res=forrCoLikelihood(c(fitE$Theta,fitE$P,fitE$rho),A,fitE$y,fitC$y[(lc-n+1):lc],fitE$opt.p,fitE$uselambda);
	}
		
	fitE$yonemu=res$yonemu	
	fitE$ssq=as.numeric(res$ssq)
	fitE$mu=res$mu
	fitE$Psi=res$Psi
	fitE$Psinv=res$Psinv
	
	fitE$nevals=nevals
	class(fitE)<- "forr"
	fit <- fitE
	fit$c <- fitC
	fit$PsidXe=res$Psi
	fit$PsinvdXe=res$Psinv
	fit$PsicXc=fitC$Psi
	fit$PsinvcXc=fitC$Psinv	
	class(fit)<- "coforr"
	fit <- forrCoModel(fit,A)
}


###################################################################################
#' Helper Function Forrester Co-Kriging
#'
#' Final calculations of the Co-Kriging fit
#'
#' @param fit co-Kriging fit of class \code{coforr}
#' @param AX differences between samples, as calculated in \code{forrCoBuilder}
#'
#' @return co-Kriging fit of class \code{coforr}, now finalized
#' @keywords internal
#' @references FORRESTER, A.I.J, SOBESTER A. & KEAN, A.J. (2007), Multi-Fidelity optimization via surrogate modelling. \emph{Proc. R. Soc. A} 463, 3251-3269. \cr
#' Forrester, Alexander I.J.; Sobester, Andras; Keane, Andy J. (2008). Engineering Design via Surrogate Modelling - A Practical Guide. John Wiley & Sons.
###################################################################################
forrCoModel <- function(fit,AX){
	Xe=fit$X
	Xc=fit$c$X
	ye=fit$y
	yc=fit$c$y
	ne=dim(Xe)[1]
	nc=dim(Xc)[1]
	#thetad=10^fit$Theta; #todo why not used?
	thetac=10^fit$c$Theta
	lambdac=10^fit$c$Lambda
	rho=fit$rho 
	y=c(yc, ye)

	
	if(fit$opt.p){ #TODO not working yet!
		AX=abs(AX)^fit$P
	}
	PsicXe=matrix(colSums(thetac*AX),ne,ne)
	PsicXe=exp(-PsicXe)+diag(lambdac,ne)
	
		
	A=matrix(0,ncol(Xe),nc*ne) #preallocate array
	for(i in 1:nc){
		A[,(1+(i-1)*ne):(ne*i)]=(Xc[i,]-t(Xe)) #TODO not a fool proof implementation #TODO improve speed?
	}		
	if(fit$opt.p){ #TODO not working yet!
		A=abs(A)^fit$P
	}else{ # p  is fixed to 2 and the array A is completely precalculated
		A=A^2
	}
	PsicXcXe=matrix(colSums(thetac*A),nc,ne,byrow=TRUE) #TODO problem somewhere over here
	PsicXcXe=exp(-PsicXcXe)+rbind(matrix(0,nc-ne,ne),diag(lambdac,ne))
	PsicXeXc=t(PsicXcXe)
	
	fit$PsicXe = PsicXe
	fit$PsicXeXc = PsicXeXc
	fit$PsicXcXe = PsicXcXe

	##########
	#build matrix covariance matrix C:
	fit$C=rbind(cbind(as.numeric(fit$c$ssq)*fit$PsicXc, 
					rho*as.numeric(fit$c$ssq)*fit$PsicXcXe),
				cbind(rho*as.numeric(fit$c$ssq)*fit$PsicXeXc, 
					rho^2*as.numeric(fit$c$ssq)*fit$PsicXe+as.numeric(fit$ssq)*fit$PsidXe))
	
	# ModelInfo.UC=chol(ModelInfo.C);
	fit$Cinv= try(solve.default(fit$C), TRUE) #TODO: this is a fix for stability, not proper math.  #Important notes: chol2inv(chol(Psi)) may be less likely to give an answer BUT may also be faster and more accurate
	if(class(fit$Cinv) == "try-error"){		#TODO: ginv may be meaningless when solve is instable...
		fit$Cinv= ginv(fit$C)  # DEPENDENCY OR SUGGEST: ginv needs MASS
	}
	fit$mu=sum(fit$Cinv%*%as.matrix(y))/sum(fit$Cinv)
	fit
}

###################################################################################
#' Calculate neg. log. lik. for Co-Kriging
#' 
#' Used to determine theta/lambda values for the Co-Kriging model in \code{\link{forrCoBuilder}}.
#'
#' @param x vector, containing log(theta) and lambda
#' @param AX 3 dimensional array, constructed by forrCoBuilder from the sample locations
#' @param Ay vector of observations at expensive sample locations
#' @param Ayc vector of observations at cheap sample locations
#'
#' @return list with elements\cr
#' \code{NegLnLike}  concentrated log-likelihood *-1 for minimising \cr
#' \code{Psi} correlation matrix\cr
#' \code{Psinv} inverse of correlation matrix \cr
#' \code{mu} \cr
#' \code{ssq}
#' @seealso \code{\link{spotPredictCoForrester}} \code{\link{forrCoBuilder}}
#' @export
#' @keywords internal
#' @references FORRESTER, A.I.J, SOBESTER A. & KEAN, A.J. (2007), Multi-Fidelity optimization via surrogate modelling. \emph{Proc. R. Soc. A} 463, 3251-3269. \cr
#' Forrester, Alexander I.J.; Sobester, Andras; Keane, Andy J. (2008). Engineering Design via Surrogate Modelling - A Practical Guide. John Wiley & Sons.
###################################################################################
forrCoLikelihood <- function(x,AX, Ay, Ayc, opt.p=FALSE,uselambda=T){
	nx<-nrow(AX)
	theta=10^x[1:nx]
	#if(opt.p){
	#	AX=abs(AX)^(x[(nx+1):(2*nx)])
	#}
	if(uselambda){
		rho=x[length(x)-1];	
		lambda=10^x[length(x)]
	}else{
		rho=x[length(x)];	
		lambda=0
	}
	if( any(c(theta)==0) ||  any(is.infinite(c(theta,lambda)))){ #unfortunately L-BFGS-B might violate bounds
		return(list(NegLnLike=1e4,Psi=NA,Psinv=NA,mu=NA,ssq=NA))
	}
	n <- dim(Ay)[1]
	#k=dim(Ay)[2] #todo why not used?
	Psi <- exp(-matrix(colSums(theta*AX),n,n))
	if(uselambda)
		Psi <- Psi+diag(lambda,n)
	# concentrated log-likelihood calculation
	LnDetPsi=as.numeric(determinant.matrix(Psi)$modulus) 
	Psinv= try(solve.default(Psi), TRUE)#Important notes: chol2inv(chol(Psi)) may be less likely to give an answer BUT may also be faster and more accurate
	if(class(Psinv) == "try-error"){		
		return(list(NegLnLike=1e4,Psi=NA,Psinv=NA,mu=NA,ssq=NA))
	}
	d=Ay - rho * Ayc #todo other models have a rho vector not a scalar value like forrester
	mu=sum(Psinv%*%d)/sum(Psinv)# note: matrix%*%onevector is signif. faster than rowSums(matrix)
	yonemu=d-mu  #%yonemu=d-one*mu
	SigmaSqr=(t(yonemu)%*%Psinv%*%yonemu)/n;
	NegLnLike=0.5*(n*log(SigmaSqr) + LnDetPsi);
	if(is.na(SigmaSqr))browser()
	list(NegLnLike=NegLnLike,Psi=Psi,Psinv=Psinv,mu=mu,yonemu=yonemu,ssq=SigmaSqr)
}

###################################################################################
#' Predict Forrester Co-Kriging Model
#' 
#' Predict new samples on a Forrester Co-Kriging model.
#'
#' @param x design matrix to be predicted 
#' @param fit fit of the Co-Kriging model (settings and parameters), as created by \code{forrCoBuilder}
#' @param pred.all if TRUE return all (RMSE and prediction, in a dataframe), else return only prediction
#'
#' @return Returned value is dependent on the setting of \code{pred.all}\cr
#' TRUE: data.frame with columns f (function values) and s (RMSE)\cr
#' FALSE: vector of function values only
#'
#' @seealso \code{\link{forrCoBuilder}}
#' @export
#' @references FORRESTER, A.I.J, SOBESTER A. & KEAN, A.J. (2007), Multi-Fidelity optimization via surrogate modelling. \emph{Proc. R. Soc. A} 463, 3251-3269. \cr
#' Forrester, Alexander I.J.; Sobester, Andras; Keane, Andy J. (2008). Engineering Design via Surrogate Modelling - A Practical Guide. John Wiley & Sons.
###################################################################################
forrCoRegPredictor <- function(x,fit,pred.all=FALSE){
	Xe=fit$X;
	Xc=fit$c$X;
	ye=fit$y;
	yc=fit$c$y;
	ne=dim(Xe)[1]
	nc=dim(Xc)[1]
	thetad=fit$dmodeltheta
	thetac=fit$c$dmodeltheta
	rho=fit$rho 
	#one=rep(1,ne+nc)
	y=c(yc, ye)
	#normalize input x
	x <- spotNormalizeMatrix2(as.matrix(x),0,1,fit$normalizexmin,fit$normalizexmax)
	Psinv=fit$Cinv #fixed: does not need to be computed, is already done in likelihood function
	mu=fit$mu
	psic=matrix(1,nrow(x),nc)
	psid=matrix(1,nrow(x),ne)
	x <- t(x) 
	nx <- ncol(x)
		
	if(nx>nc){ #loop over whatever is smaller: observed samples or samples to be predicted
		for (i in 1:nc)
			psic[,i]=colSums(thetac*((Xc[i,]-x)^2))		
	}else{			
		psic <- NULL
		for (i in 1:nx)
			psic=rbind(psic,colSums(thetac*(t(Xc)-x[,i])^2))
	}			
	psic = rho*fit$c$ssq*exp(-psic)
	
	if(nx>ne){ #loop over whatever is smaller: observed samples or samples to be predicted
		for (i in 1:ne){
			diffX <- (Xe[i,]-x)^2
			psid[,i]=rho^2*fit$c$ssq*exp(-colSums(thetac*diffX))+fit$ssq*exp(-colSums(thetad*diffX))
		}
	}else{			
		psid <- NULL
		for (i in 1:nx){
			diffX <- (t(Xe)-x[,i])^2
			psid=rbind(psid,rho^2*fit$c$ssq*exp(-colSums(thetac*diffX))+fit$ssq*exp(-colSums(thetad*diffX)))
		}	
	}	

	psi <- cbind(psic,psid)
	f=as.numeric(psi%*%(Psinv%*%(y-mu)))+mu #vectorised
	##########################################################################
	if (pred.all){
		lambda=fit$dmodellambda; 
		SSqr= rho^2*fit$c$ssq+fit$ssq+lambda-diag(psi%*%(Psinv%*%t(psi))) #vectorised
		s=sqrt(abs(SSqr));
	}
	if(!pred.all){list(f=f)}else{data.frame(f=f,s=as.numeric(s))}
}
###################################################################################
#' Meta Model Interface: DACE Kriging
#' 
#' This Kriging meta model is based on DACE (Design and Analysis of Computer Experiments).
#' It allows to choose different regression and correlation models. If multiple response 
#' variables are present, a DACE model for each will be created for the purpose of multi objective optimization.
#'
#' @param rawB unmerged data 
#' @param mergedB merged data (aggregation of repeated design points)
#' @param design new design points which should be predicted
#' @param spotConfig global list of all options, needed to provide data for calling functions. This also contains a list, with settings for Forrester:\cr
#'	\code{spotConfig$seq.dace.budget} Budget for MLE of parameters, default is 100. The value will be multiplied with the length of the model parameter vector to be optimized. (which means 100*i evaluations, where i depends on the problem dimension as well as the correlation function and whether a nugget value is to be determined)\cr
#'	\code{spotConfig$seq.dace.tol} Tolerance stopping criterion for MLE, default is 1e-6\cr
#'	\code{spotConfig$seq.dace.regr} Regression function to be used: \code{\link{regpoly0}} (default), \code{\link{regpoly1}}, \code{\link{regpoly2}}\cr
#'	\code{spotConfig$seq.dace.corr} Correlation function to be used: \code{\link{corrnoisykriging}} (default), \code{\link{corrkriging}}, \code{\link{corrnoisygauss}}, \code{\link{corrgauss}}, \code{\link{correxp}}, \code{\link{correxpg}}, \code{\link{corrlin}}, \code{\link{corrcubic}},\code{\link{corrspherical}},\code{\link{corrspline}}\cr
#'	\code{spotConfig$seq.dace.nugget} Value for nugget. Default is -1, which means the nugget will be optimized during MLE. Else it can be fixed in a range between 0 and 1.
#' @param fit if an existing model fit is supplied, the model will not be build based on 
#'				data, but only evaluated with the model fit (on the design data). To build the model, 
#'				this parameter has to be NULL. If it is not NULL the parameters mergedB and rawB will not be 
#'				used at all in the function.
#'
#' @return returns the list \code{spotConfig} with two new entries:\cr
#' 	spotConfig$seq.modelFit fit of the model used with dacePredictor \cr
#'	spotConfig$seq.largeDesignY the y values of the design, evaluated with the fit
#'
#' @seealso \code{\link{dacePredictor}} \code{\link{daceBuilder}}
#'
#' @author The authors of the original DACE Matlab toolbox \url{http://www2.imm.dtu.dk/~hbni/dace/} 
#' are Hans Bruun Nielsen \email{hbn@@imm.dtu.dk}, Soren Nymand Lophaven and Jacob Sondergaard. \cr
#' Extension of the Matlab code by Tobias Wagner \email{wagner@@isf.de}. \cr 
#' Porting and adaptation to R and further extensions by Martin Zaefferer \email{martin.zaefferer@@fh-koeln.de}.
#'
#' @references S.~Lophaven, H.~Nielsen, and J.~Sondergaard.
#' {DACE---A Matlab Kriging Toolbox}.
#' Technical Report IMM-REP-2002-12, Informatics and Mathematical
#' Modelling, Technical University of Denmark, Copenhagen, Denmark, 2002.
#'
#' @export
###################################################################################
spotPredictDace <- function(rawB,mergedB,design,spotConfig,fit=NULL){	#TODO: Implement smoothing
	design <- spotInitializePredictor(design,"data.frame",spotConfig$alg.roi,NULL,"spotPredictDace",spotConfig$io.verbosity)
	
	########################################################
	# BUILD
	########################################################	
	if(is.null(fit)){
		
		########################################################
		# Parameter defaults
		########################################################	
		if(is.null(spotConfig$seq.dace.budget))spotConfig$seq.dace.budget=100; #budget for lik. estimation
		if(is.null(spotConfig$seq.dace.tol))spotConfig$seq.dace.tol=1e-6; #tolerance for lik. estimation
		if(is.null(spotConfig$seq.dace.regr))spotConfig$seq.dace.regr=regpoly0; #regression function
		if(is.null(spotConfig$seq.dace.nugget))spotConfig$seq.dace.nugget=-1; #-1 means nugget is optimized, else it should be in inteval between 0 and 1
		if(is.null(spotConfig$seq.dace.corr))spotConfig$seq.dace.corr=corrnoisykriging; #-1 means nugget is optimized, else it should be in inteval between 0 and 1
		if(is.null(spotConfig$seq.dace.algtheta))spotConfig$seq.dace.algtheta="optim-L-BFGS-B"; #optimization algorithm to be used
		
		xNames <- row.names(spotConfig$alg.roi);
		yNames <-  spotConfig$alg.resultColumn
		x <- as.matrix(rawB[xNames])
		if(length(yNames)==1){
			y <- rawB[[yNames]]
			fit<-daceBuilder(x, y, tol=spotConfig$seq.dace.tol, budget=spotConfig$seq.dace.budget, regr=spotConfig$seq.dace.regr, corr=spotConfig$seq.dace.corr, nugget= spotConfig$seq.dace.nugget,algtheta=spotConfig$seq.dace.algtheta)
		}
		else{#Distinction for multi criteria spot 
			y <- rawB[yNames]
			fit=list()			
			for (i in 1:length(yNames)){
				fit[[i]]<-daceBuilder(x, y[,i], tol=spotConfig$seq.dace.tol, budget=spotConfig$seq.dace.budget, regr=spotConfig$seq.dace.regr, corr=spotConfig$seq.dace.corr, nugget= spotConfig$seq.dace.nugget,algtheta=spotConfig$seq.dace.algtheta)
			}			
		}
	}
	########################################################
	# PREDICT
	########################################################	
	if(!is.null(design)){ 	
		MSE=spotConfig$seq.model.variance
		nmodel <- length(spotConfig$alg.resultColumn)
		if(nmodel>1){ #do multi criteria prediction
			resy=matrix(0,nrow(design),nmodel)
			resvar=matrix(NA,nrow(design),nmodel)
			y=list()
			for (i in 1:length(fit)){ #predict			
				res= dacePredictor(as.matrix(design), fit[[i]],GRAD=FALSE,MSE=MSE,GRADMSE=FALSE)
				resy[,i]= res$f
				if(MSE)resvar[,i]= res$s
			}
			y=spotConfig$alg.currentResult[spotConfig$alg.resultColumn] #TODO better from mergedB?
			if(is.function(spotConfig$seq.infill)){# do EI 
				resy= spotConfig$seq.infill(resy,resvar,y,spotConfig$mco.refPoint)
			}
		}else{ #do single criteria prediction
			res= dacePredictor(as.matrix(design), fit,GRAD=FALSE,MSE=MSE,GRADMSE=FALSE)
			resy=res$f
			resvar=matrix(NA,nrow(design),1)
			if(MSE)resvar=res$s
			if(is.function(spotConfig$seq.infill)){ # do EI			
				resy= spotConfig$seq.infill(resy,resvar,min(spotConfig$alg.currentResult[spotConfig$alg.resultColumn])) #TODO better from mergedB?
			}
		}
	}else{
		resy <- NULL
		resvar <- NULL
	}	
	########################################################
	# Output
	########################################################	
	spotWriteLines(spotConfig$io.verbosity,3,"spotPredictDace finished");	
	spotConfig$seq.largeDesignY=as.data.frame(resy)
	spotConfig$seq.largeDesignVar=as.data.frame(resvar)	
	spotConfig$seq.modelFit<-fit;
	spotConfig
}

###################################################################################
#' DACE predictor
#' 
#' Predicts y(x) for a given DACE model (i.e. as created by \code{\link{daceBuilder}}).
#' The user can choose whether to predict only mean or if he is also interested in gradient, mean squared error MSE, or the MSE gradient.
#'
#' @param x \code{mx} candidate points of dimension \code{n} to be predicted. That is, a matrix with \code{mx} rows and \code{n} columns.
#' @param fit the model fit, as returned by \code{\link{daceBuilder}}
#' @param GRAD specify whether gradient of response should be computed (default is FALSE, which means no). Even if GRAD is TRUE, the gradient will only be computed in case of a single design point \code{mx==1}.
#' @param MSE specify whether estimated MSE of response should be computed (default is FALSE, which means no)
#' @param GRADMSE specify whether gradient of MSE should be computed (default is FALSE, which means no).  Even if GRADMSE is TRUE, the gradient will only be computed in case of a single design point \code{mx==1}.
#'
#' @return returns a list with the following elements: 
#' 			\item{\code{f}}{ Predicted response \code{y} at design points \code{x} (always)}
#' 			\item{\code{df}}{ Gradient of response \code{y} at design points \code{x}  (only if: \code{GRAD==TRUE} and \code{mx==1})}
#' 			\item{\code{s}}{ Estimated MSE (only if: \code{MSE==TRUE})}
#' 			\item{\code{ds}}{ Gradient of MSE (only if: \code{GRADMSE==TRUE} and \code{mx==1})}
#'
#' @seealso \code{\link{spotPredictDace}} \code{\link{daceBuilder}}
#'
#' @author The authors of the original DACE Matlab toolbox \url{http://www2.imm.dtu.dk/~hbni/dace/} 
#' are Hans Bruun Nielsen \email{hbn@@imm.dtu.dk}, Soren Nymand Lophaven and Jacob Sondergaard. \cr
#' Additional code for generalization to different models by Tobias Wagner \email{wagner@@isf.de}. \cr 
#' Porting and adaptation to R and further extensions by Martin Zaefferer \email{martin.zaefferer@@fh-koeln.de}.
#'
#' @references S.~Lophaven, H.~Nielsen, and J.~Sondergaard.
#' {DACE---A Matlab Kriging Toolbox}.
#' Technical Report IMM-REP-2002-12, Informatics and Mathematical
#' Modelling, Technical University of Denmark, Copenhagen, Denmark, 2002.
#'
#' @examples
#' ## Create design points
#' x = cbind(runif(20)*15-5,runif(20)*15)
#' ## Compute observations at design points (for Branin function)
#' y = apply(x,1,spotBraninFunction)
#' ## Create model
#' fit = daceBuilder(x,y)
#' ## Create candidate design points
#' xx = cbind(runif(20)*15-5,runif(20)*15)
#' ## Predict candidates
#' y1 = dacePredictor(xx,fit)$f
#' ## Plot residuals
#' plot(y1 - apply(xx,1,spotBraninFunction))
#' ## Plot model (in comments, due to time consumption)
#' #fn <- function(x){dacePredictor(as.matrix(x),fit)$f}
#' #spotSurf3d(fn,c(-5,0),c(10,15))
#' ## Plot real function
#' #spotSurf3d(function(x){apply(x,1,spotBraninFunction)},c(-5,0),c(10,15))
#'
#' @export
###################################################################################
dacePredictor = function (x, fit,GRAD=FALSE,MSE=FALSE,GRADMSE=FALSE){
	dmodel=fit$model
	if(GRADMSE && !MSE)stop("To evaluate gradient of mse, both GRADMSE and MSE must be true.")
	mse = NaN;   grad = NaN;  dmse = NaN;  #% Default return values
	sx = dim(x);            #% number of trial sites and their dimension
	if(is.null(sx))sx=c(1,length(x));
	if (sx[1] > 5000){
		y =rep(0,sx[1]);
		mse = rep(0,sx[1]);
		for(i in 1:floor(sx[1]/5000)){
			res=dacePredictor(x[((i-1)*5000+1):(i*5000),], fit,FALSE,MSE,FALSE)
			y[((i-1)*5000+1):(i*5000)]= res$f
			if(MSE)	mse[((i-1)*5000+1):(i*5000)]= res$s;
		}
		if(!(i*5000+1)>sx[1]){
			res=dacePredictor(x[(i*5000+1):sx[1],], fit)
			y[(i*5000+1):sx[1]]= res$f
		}
		if(MSE){
			if(!(i*5000+1)>sx[1]){mse[(i*5000+1):sx[1]]= res$s}
			return(data.frame(f=as.vector(y),s=mse))
		}else{return(data.frame(f=as.vector(y)))}
	}

	if(any(is.nan(dmodel$beta)))	stop('DMODEL has no beta value');	

	mn = dim(dmodel$S);  #% number of design sites and number of dimensions
	if(is.null(mn)){
		m=1
		n=length(dmodel$S)
	}else{
		m=mn[1];
		n=mn[2];
	}
	
	
	
	#if( min(sx) == 1 && n > 1){ #% Single trial point. 
	#	nx = max(sx);
	#	if  (nx == n){
	#		mx = 1;  
	#		x = t(x); #TODO?
	#	}
	#}else{
	#	mx = sx[1];  nx = sx[2];
	#}
	#### better
	mx=sx[1]
	nx=sx[2]
	if(nx != n) stop(paste('Dimension of trial sites should be',n))


	#% Normalize trial sites
	#x = (x - t(repmat(dmodel$Ssc[1,],1,mx))) /  t(repmat(dmodel$Ssc[2,],1,mx));   
	x = (x - dmodel$Ssc[rep(1,mx),] )/   dmodel$Ssc[rep(2,mx),];	#MZ: fix for faster processing
	#q = size(dmodel.Ysc,2);  #% number of response functions
	q = 1;   #TODO adapt dace for multiple responses? other functions need to be tested for this first
	#y = zeros(mx,q);         #% initialize result
	y = rep(0,mx);         	
	
	if (mx == 1){  #% one site only
		#dx = repmat(x, m, 1) - dmodel$S;  #% distances to design sites
		dx = matrix(x,ncol=length(x),nrow=m,byrow=TRUE) - dmodel$S
		res = dmodel$regr(x,TRUE);
		f = res$f 
		df=res$df
		if(GRAD || GRADMSE){
			res = dmodel$corr(dmodel$theta, dx, "all");
			r = res$r 
			dr=res$dr
		}else{
			res = dmodel$corr(dmodel$theta, dx, "r");
			r = res$r 		
		}
		#% Scaled predictor
		sy = f %*% dmodel$beta + t(dmodel$gamma%*%r);
		#% Predictor
		y = t(dmodel$Ysc[1] + dmodel$Ysc[2] * sy);  #todo this needs adaption if q not 1
		results=data.frame(f=as.vector(y))
		if (GRAD){               #% gradient/Jacobian wanted
			#% Scaled Jacobian
			dy = df * dmodel$beta + dmodel$gamma %*% dr;  #todo first term matmult?
			#% Unscaled Jacobian
			grad = dy * repmat(dmodel$Ysc[2], 1, nx) / t(repmat(dmodel$Ssc[2], 1, q)); #todo this needs adaption if q not 1, or switch everything to fixed q=1 ?   #TODO remove repmat
			results$df=grad
		}
		if(MSE){  #% MSE wanted				
			rt = spotHelpBslash(dmodel$C,r);
			u = t(dmodel$Ft) %*% rt - t(f); 
			v = spotHelpBslash(dmodel$G,u);
			#mse = t(repmat(dmodel$sigma2,1,mx)) * repmat((1 + sum(v^2) - sum(rt^2)),1,q);	
			mse = dmodel$sigma2 * (1 + sum(v^2) - sum(rt^2));	 #TODO q was completely removed here...
			results$s=abs(mse) #TODO sometimes this seems to become a very small negative number... numerical problem? rather cut to zero instead of abs? see also below.
			if (GRADMSE){  #% gradient/Jacobian of MSE wanted
				#% Scaled gradient as a row vector
				Gv = spotHelpBslash(t(dmodel$G),v);
				g = t((dmodel$Ft %*% Gv - rt)) %*% spotHelpBslash(dmodel$C,dr) - t(df * Gv); # last term matmult??
				#% Unscaled Jacobian
				dmse = repmat(2 * dmodel$sigma2,1,nx) * (repmat(g / dmodel$Ssc[2,],1,q)); #TODO remove repmat
				results$ds=dmse
			}
		}		
	}else{ # % several trial sites
		#% Get distances to design sites
		dx = matrix(0,mx*m,n);  kk = 1:m;
		for (k in 1:mx){
			dx[kk,] = x[rep(k,m),] - dmodel$S;
			kk = kk + m;
		}
		#% Get regression function and correlation
		f = dmodel$regr(x)$f 
		r = dmodel$corr(dmodel$theta, dx, "r")$r 
		r = matrix(r,nrow=m);   

		#% Scaled predictor
		sy = f %*% dmodel$beta + t(dmodel$gamma %*% r);
		#% Predictor
		#y = t(repmat(dmodel$Ysc[1],1,mx)) + t(repmat(dmodel$Ysc[2],1,mx)) * sy;		 #todo this needs adaption if q not 1, or switch everything to fixed q=1 ?
		y = dmodel$Ysc[1] + dmodel$Ysc[2] * sy;		 #todo this needs adaption if q not 1, or switch everything to fixed q=1 ?
		results=data.frame(f=as.vector(y))
		if (MSE) {  #% MSE wanted
			rt = spotHelpBslash(dmodel$C,r);
			u = as.matrix(spotHelpBslash(dmodel$G,(t(dmodel$Ft) %*% rt - t(f))));
			#mse = t(repmat(dmodel$sigma2,1,mx)) * repmat((1 + rowSums(u^2) - colSums(rt^2)),1,q);
			mse = dmodel$sigma2 * (1 + colSums(u^2) - colSums(rt^2));    #q removed for more efficient code
			if(GRAD || GRADMSE)	warning('Only  y  and  mse  are computed for multiple design sites');
			results$s=abs(mse)
		}
	} #% of several sites
	results
}
###################################################################################
#' Print Function DACE Kriging
#'
#' Print information about a DACE Kriging fit, as produced by \code{\link{daceBuilder}}.
#'
#' @rdname print
#' @method print dace
# @S3method print dace
#' @param x	fit returned by \code{\link{forrBuilder}}.
#' @param ... additional parameters	
#' @export
#' @keywords internal
###################################################################################
print.dace <- function(x,...){
	cat("------------------------\n")
	cat("Dace Kriging model.\n")
	cat("------------------------\n")
	cat("Estimated activity parameters (theta) sorted \n")
	cat("from most to least important variable \n")
	cat(paste("x",order(x$theta,decreasing=TRUE),sep="",collaps=" "))
	cat("\n")	
	cat(sort(x$theta,decreasing=TRUE))	
	cat("\n \n")
	cat("exponent(s) p:\n")
	if(!is.null(x$p))
		cat(x$p)
	else
		cat("none (no exponents in the employed correlation function)")
	cat("\n \n")
	cat("Estimated regularization constant (or nugget) lambda:\n")	
	if(!is.null(x$p))
		cat(x$lambda)
	else
		cat("none (no nugget in the employed correlation function)")
	cat("\n \n")
	cat("Number of Likelihood evaluations during MLE:\n")	
	cat(x$nevals)
	cat("\n")	
	cat("------------------------\n")
}

###################################################################################
#' Build DACE model
#'
#' This Kriging meta model is based on DACE (Design and Analysis of Computer Experiments).
#' It allows to choose different regression and correlation models. The optimization of model parameters
#' is by default done with a bounded simplex method from the \code{nloptr} package.
#'
#' @param parameters known design points. That is, a matrix with \code{n} rows (for each point) and \code{dim} columns (for each dimension).
#' @param objectives vector of observations at known design points of length \code{n}.
#' @param startTheta vector which contains initial guess for model parameters theta. Initial guess will be set depending on correlation function if this vector is missing.
#' @param tol Tolerance stopping criterion for the simplex MLE. Default is \code{1e-6}.
#' @param budget Number of Likelihood evaluations, as a stopping criterion for the simplex MLE. Default is \code{Inf}. The value will be multiplied with the length of the model parameter vector to be optimized.
#' @param regr Regression function to be used: \code{\link{regpoly0}} (default), \code{\link{regpoly1}}, \code{\link{regpoly2}}. Can be a custom user function.
#' @param corr Correlation function to be used: \code{\link{corrnoisykriging}} (default), \code{\link{corrkriging}}, \code{\link{corrnoisygauss}}, \code{\link{corrgauss}}, \code{\link{correxp}}, \code{\link{correxpg}}, \code{\link{corrlin}}, \code{\link{corrcubic}},\code{\link{corrspherical}},\code{\link{corrspline}}. Can also be user supplied (if in the right form).
#' @param nugget Value for nugget. Default is -1, which means the nugget will be optimized during MLE. Else it can be fixed in a range between 0 and 1.
#' @param algtheta algorithm used to find theta, default is \code{"optim-L-BFGS-B"}. Else, any from the list of possible \code{method} values in \code{\link{spotOptimizationInterface}} can be chosen.
#'
#' @return returns a list with the following elements:
#' 			\item{\code{model}}{ Again a list, containing model parameters}
#' 			\item{\code{like}}{ Likelihood value}
#' 			\item{\code{theta}}{activity parameters theta (vector)}
#' 			\item{\code{p}}{exponents p (vector)}
#' 			\item{\code{lambda}}{nugget value (numeric)}
#' 			\item{\code{nevals}}{ Number of iterations during MLE}
#'
#' @seealso \code{\link{spotPredictDace}} \code{\link{dacePredictor}}
#'
#' @author The authors of the original DACE Matlab toolbox \url{http://www2.imm.dtu.dk/~hbni/dace/} 
#' are Hans Bruun Nielsen \email{hbn@@imm.dtu.dk}, Soren Nymand Lophaven and Jacob Sondergaard. \cr
#' Extension of the Matlab code by Tobias Wagner \email{wagner@@isf.de}. \cr 
#' Porting and adaptation to R and further extensions by Martin Zaefferer \email{martin.zaefferer@@fh-koeln.de}.
#'
#' @references S.~Lophaven, H.~Nielsen, and J.~Sondergaard.
#' {DACE---A Matlab Kriging Toolbox}.
#' Technical Report IMM-REP-2002-12, Informatics and Mathematical
#' Modelling, Technical University of Denmark, Copenhagen, Denmark, 2002.
#'
#' @examples
#' ## Create design points
#' x = cbind(runif(20)*15-5,runif(20)*15)
#' ## Compute observations at design points (for Branin function)
#' y = apply(x,1,spotBraninFunction)
#' ## Create model with default settings
#' fit = daceBuilder(x,y)
#' ## Print model parameters
#' print(fit)
#' ## Create with different regression and correlation functions
#' fit = daceBuilder(x,y,regr=regpoly2,corr=corrspline)
#' ## Print model parameters
#' print(fit)
#'
#' @export
###################################################################################
daceBuilder <- function(parameters, objectives, startTheta, tol=1e-6, budget=100 ,regr=regpoly0, corr= corrnoisykriging ,nugget = -1, algtheta="optim-L-BFGS-B"){ #nugget -1 means that the nugget will be optimized in lme
	size= dim(parameters);
	n= size[1]
	m= size[2]
	res <- daceStartParameters(n,m,nugget,corr)
	lb <- res$lb
	ub <- res$ub 
	if(missing(startTheta)){
		startTheta = res$theta; 
	}	
	
	q = length(startTheta)*budget;
	#q = min(q,budget)
	
	small<-startTheta<lb  #force starting guess into bounds
	big<-startTheta>ub
	startTheta[small]<-lb[small]
	startTheta[big]<-ub[big]	
	
	para <- dacePrepareFit(parameters, objectives, nugget, regr, corr) 
	opts=list(fevals=q, reltol=tol, restarts=FALSE)	#TODO false? TODO make restarts a parameter?
	res <- spotOptimizationInterface(par=startTheta,fn=daceLikelihood,gr=NULL,lower=lb,upper=ub,method=algtheta,
						control=opts,para=para,nugget=nugget)				
			
	bestTheta = res$par		

	ftheta <- daceFixTheta(m,bestTheta,nugget,corr)
	model = daceGetFit(ftheta$thetaConv,para)
	like =  res$value
	nEval = as.numeric(res$counts[[1]])
	fit <- list(model=model,like=like,theta=ftheta$theta,lambda=ftheta$lambda,p=ftheta$p,nevals=nEval)
	class(fit)<- "dace"
	fit
}


###################################################################################
#' Start parameter setup
#'
#' This function returns a starting guess for theta, as well as suitable lower and upper bounds.
#' The result depends on dimensionality of the problem, number of design points, the nugget value and the choice of correlation function.
#'
#' @param n number of known design points
#' @param m dimension (length) of each point
#' @param nugget Value for nugget. Default is -1, which means the nugget will be optimized during MLE. In that case, a lower limit of 0.5 and an upper limit of 1 as well as a starting value of 0.999 will added to the three output vectors (theta, lower and upper bounds). This is only relevant for correlation functions that use a nugget (\code{\link{corrnoisygauss}},\code{\link{corrnoisykriging}} )
#' @param corr The choice of correlation function (which defines the length and values of theta and bounds): \code{\link{corrnoisykriging}} (default), \code{\link{corrkriging}}, \code{\link{corrnoisygauss}}, \code{\link{corrgauss}}, \code{\link{correxp}}, \code{\link{correxpg}}, \code{\link{corrlin}}, \code{\link{corrcubic}},\code{\link{corrspherical}},\code{\link{corrspline}}. Can also be user supplied (if in the right form).
#'
#' @return returns a list with the following elements: \cr
#' 			\code{theta} Starting point for the internal parameter estimation\cr
#' 			\code{lb} lower bound\cr
#' 			\code{ub} upper bound
#'
#' @seealso \code{\link{daceBuilder}}
#' @keywords internal
# @export
###################################################################################
daceStartParameters  <- function(n,m,nugget,corr){
	if(identical(corr,corrnoisykriging)||identical(corr,corrkriging)){
		ones = matrix(1,1,m)
		lb = c(ones*(-12), ones*0.01);   #TODO wouldnt rep(-12,m) be faster? or would that yield different class ? (vec, mat)
		ub = c(ones*10, ones*2);
		theta = c(ones*(n/(100*m)), 1.9*ones)
	#}else if(identical(corr,corrnoisygauss)||identical(corr,corrspline)||identical(corr,corrspherical)||identical(corr,corrgauss)||identical(corr,corrlin)||identical(corr,corrcubic)||identical(corr,correxp)){ # like corrnoisykriging, but exponents not optimized, only activity parameters
	}else if(identical(corr,correxpg)){
		ones = matrix(1,1,m)
		lb = c(ones*(-12),-6);   
		ub = c(ones*10,2);
		theta = c(ones*(n/(100*m)),-1)
	}else {	
		ones = matrix(1,1,m)
		lb = ones*(-12);   
		ub = ones*10;
		theta = ones*(n/(100*m))
	}	
	
	if(nugget==-1 && (identical(corr,corrnoisygauss) || identical(corr,corrnoisykriging))){#nugget is optimized, else fixed
		lb = c(lb, 0.5)
		ub = c(ub, 1)
		theta = c(theta,0.999)
	}
	
	list(theta=theta,lb=lb,ub=ub)
}

###################################################################################
#' Fix model parameters
#'
#' This function fixes theta (model parameter vector) after optimization. That is, if a value was optimized on a logarithmic scale: \code{theta=10^theta}.
#' The result depends on dimensionality of the problem, the nugget value and the choice of correlation function.
#'
#' @param m dimension (length) of each point
#' @param bestTheta the theta value found during MLE
#' @param nugget Value for nugget. Default is -1, which means the nugget was be optimized during MLE. In this case, it is part of bestTheta. Otherwise, \code{nugget} is appended to the theta vector. This is only relevant for correlation functions that use a nugget (\code{\link{corrnoisygauss}},\code{\link{corrnoisykriging}} )
#' @param corr The choice of correlation function (which defines the length and values of theta and bounds): \code{\link{corrnoisykriging}} (default), \code{\link{corrkriging}}, \code{\link{corrnoisygauss}}, \code{\link{corrgauss}}, \code{\link{correxp}}, \code{\link{correxpg}}, \code{\link{corrlin}}, \code{\link{corrcubic}},\code{\link{corrspherical}},\code{\link{corrspline}}. Can also be user supplied (if in the right form).
#'
#' @return returns a list:
#'	\code{thetaConv} the fixed theta vector (all model parameters)
#'	\code{theta} vector of activity parameters theta
#'	\code{p} vector of exponents p(NULL if not used in correlation function)
#'	\code{lambda} nugget value lambda (NULL if not used in correlation function)
#'
#' @seealso \code{\link{daceBuilder}} 
#' @keywords internal
# @export
###################################################################################
daceFixTheta  <- function(m,bestTheta,nugget,corr){
	lambda <- NULL
	p <- NULL
	if(identical(corr,corrnoisykriging)){
		if(nugget > 0 && nugget <= 1){
			thetaConv <- c(10^bestTheta[1:m], bestTheta[(m+1):(2*m)], nugget)
			theta <- thetaConv[1:m]
			p <- thetaConv[(m+1):(2*m)]
			lambda <- nugget
		}else if(nugget==-1){
			thetaConv <- c(10^bestTheta[1:m], bestTheta[(m+1):(2*m+1)])
			theta <- thetaConv[1:m]
			p <- thetaConv[(m+1):(2*m)]
			lambda <- thetaConv[2*m+1]			
		}
	}else if(identical(corr,corrkriging)){
		thetaConv <- c(10^bestTheta[1:m], bestTheta[(m+1):(2*m)])
		theta <- thetaConv[1:m]
	}else if(identical(corr,corrnoisygauss)){ # like corrnoisykriging, but exponents not optimized, only activity parameters
		if(nugget > 0 && nugget <= 1){
			thetaConv <- c(10^bestTheta[1:m], nugget)
			theta <- thetaConv[1:m]
			p <- 2
			lambda <- nugget			
		}else if(nugget==-1){
			thetaConv <- c(10^bestTheta[1:m], bestTheta[m+1])
			theta <- thetaConv[1:m]
			p <- 2
			lambda <- thetaConv[m+1]			
		}
	#}else if(identical(corr,correxpg)||identical(corr,corrgauss)||identical(corr,corrspline)||identical(corr,corrspherical)||identical(corr,corrlin)||identical(corr,corrcubic)||identical(corr,correxp)){ # like corrnoisykriging, but exponents not optimized, only activity parameters
	}else{
		thetaConv <- 10^bestTheta
		theta <- thetaConv
	}
	list(thetaConv=thetaConv,theta=theta,lambda=lambda,p=p)
}

###################################################################################
#' Wrapper for Maximum Likelihood Estimation
#'
#' Returns the maximum likelihood for the model parameter optimization.
#'
#' @param theta model parameter vector to be evaluated
#' @param para model option list, as created with \code{\link{dacePrepareFit}}
#' @param nugget Value for nugget. Default is -1, which means the nugget was optimized during MLE. 
#'
#' @return the likelihood value as calculated by \code{\link{daceEvalFit}}
#'
#' @seealso \code{\link{daceBuilder}} \code{\link{daceEvalFit}} 
#' @keywords internal
# @export
###################################################################################
daceLikelihood <- function (theta, para, nugget){
	n <- para$n
	thetaConv <- daceFixTheta(n,theta,nugget,para$corr)$thetaConv
	res <- daceEvalFit(thetaConv,para)		
	min(res[length(thetaConv)+1])  #TODO: warum hier min?
}


###################################################################################
#' Prepare DACE fit
#'
#' Prepares a list with relevant model options and settings based on user choice and problem setup.
#'
#' @param S known design points. That is, a matrix with \code{n} rows (for each point) and \code{dim} columns (for each dimension).
#' @param Y vector of observations at known design points of length \code{n}.
#' @param nugget Value for nugget. Default is -1, which means the nugget will be optimized during MLE. 
#' @param regr Regression function to be used: \code{\link{regpoly0}} (default), \code{\link{regpoly1}}, \code{\link{regpoly2}}. Can be a custom user function.
#' @param corr Correlation function to be used: \code{\link{corrnoisykriging}} (default), \code{\link{corrkriging}}, \code{\link{corrnoisygauss}}, \code{\link{corrgauss}}, \code{\link{correxp}}, \code{\link{correxpg}}, \code{\link{corrlin}}, \code{\link{corrcubic}},\code{\link{corrspherical}},\code{\link{corrspline}}. Can also be user supplied (if in the right form).
#'
#' @return a list with several model or problem specific settings and parameters
#'
#' @author The authors of the original DACE Matlab code are Hans Bruun Nielsen \email{hbn@@imm.dtu.dk}, Soren Nymand Lophaven and Jacob Sondergaard. \cr
#' Extension of the Matlab by Tobias Wagner \email{wagner@@isf.de}. \cr 
#' Porting and adaptation to R and further extensions by Martin Zaefferer \email{martin.zaefferer@@fh-koeln.de}.
#'
#' @seealso \code{\link{daceBuilder}} 
#' @keywords internal
# @export
###################################################################################
dacePrepareFit <- function (S, Y, nugget, regr=regpoly0, corr=corrnoisykriging){	#this was cut from the dacefit function. dacefit was called repeatedly, doing this same thing again and again. now only done once or twice.
	m <- dim(S)[1]
	n <- dim(S)[2]
	sY <- length(Y)
	lY=sY[1]
	if (m != lY){stop("S (design) and Y (observations) must have the same number of rows")}
	#% Normalize data
	mS = colMeans(S)
	sS = apply(S,2,sd)
	mY = mean(Y)
	sY = sd(Y)

	#% 02.08.27: Check for 'missing dimension'
	sS[sS==0] = 1
	sY[sY==0] = 1
	
	S = (S - matrix(mS,ncol=length(mS),nrow=m,byrow=TRUE)) / matrix(sS,ncol=length(sS),nrow=m,byrow=TRUE)

	Y = as.matrix((Y - mY) / sY)
	#% Calculate distances D between points
	mzmax = m*(m-1) / 2        #% number of non-zero distances
	ij =  matrix(0,mzmax,2)       #% initialize matrix with indices
	D = matrix(0,mzmax,n)        #% initialize matrix with distances
	ll = 0
	for (k in 1:(m-1)){
		ll = tail(ll,1) + (1 : (m-k))
		ij[ll,] = c(rep(k, m-k),(k+1):m) 	
		D[ll,] =  S[rep(k,m-k),] - S[(k+1):m,] #% differences between points  		#MZ: faster than repmat. do similar fixes elsewhere TODO
	}	
	if( !(identical(corr,corrnoisykriging) || identical(corr,corrnoisygauss)) &  (min(sum(abs(D),2))==0)){
		stop('Multiple design sites are not allowed')
	}
	#% Regression matrix
	F = regr(S)$f 
	dims= dim(F)
	if(is.null(dims)){
		mF = length(F)
		p=1
		if  (mF != m) stop('number of rows in  F  and  S  do not match')
	}else{
		mF = dims[1]
		p = dims[2]
		if  (mF != m) stop('number of rows in  F  and  S  do not match')
		if  (p > mF) stop('least squares problem is underdetermined')
	}

	#% parameters for objective function
	list(n=n, corr=corr, regr=regr, y=Y, F=F, D=D, ij=ij, sY=sY, scS=sS, Ysc=c(mY, sY),Ssc=rbind(mS, sS),S=S)
}


###################################################################################
#' Evaluate DACE fit
#'
#' Evaluate the fit of a certain set of model parameters (\code{theta}).
#'
#' @param theta model parameters to be evaluated
#' @param para model option list, as created with \code{\link{dacePrepareFit}}
#'
#' @return performance vector, first elements are theta, last element is likelihood.
#'
#' @seealso \code{\link{daceBuilder}} \code{\link{daceLikelihood}} \code{\link{daceGetFit}} 
#' @keywords internal
# @export
###################################################################################
daceEvalFit <- function (theta, para){	
	if(any(theta <= 0)){ 
		stop('theta for dace model in spotPredictDace must be strictly positive')
	}

	f = daceObjfunc(theta, para, "f")$f
	perf = c(theta, f, 1)
	#if(is.infinite(f)) warning('Bad point.  Try increasing theta');

	perf
}  

###################################################################################
#' Get DACE fit
#'
#' Evaluate the fit of a certain set of model parameters (\code{theta}), and
#' get all relevant variables of the model.
#'
#' @param theta model parameters to be evaluated
#' @param para model option list, as created with \code{\link{dacePrepareFit}}
#'
#' @return list of model variables, with the following elements: \cr
#' 		\code{regr} regression function used \cr
#' 		\code{corr} correlation function used \cr
#' 		\code{theta} model parameters\cr
#' 		\code{beta} generalized least squares estimate\cr
#' 		\code{gamma} correlation factors\cr
#' 		\code{sigma2} maximum Likelihood estimate of the process variance\cr
#' 		\code{S} scaled design points\cr
#' 		\code{Ssc} scaling factors for design arguments\cr
#' 		\code{Y} scaled observations\cr
#' 		\code{Ysc} scaling factors for design ordinates\cr
#' 		\code{C}  Cholesky factor of correlation matrix\cr
#' 		\code{Ft} Decorrelated regression matrix\cr
#' 		\code{G} From QR factorization: Ft = Q*t(G)\cr
#'
#' @author The authors of the original DACE Matlab code are Hans Bruun Nielsen \email{hbn@@imm.dtu.dk}, Soren Nymand Lophaven and Jacob Sondergaard. \cr
#' Extension of the Matlab code by Tobias Wagner \email{wagner@@isf.de}. \cr 
#' Porting and adaptation to R and further extensions by Martin Zaefferer \email{martin.zaefferer@@fh-koeln.de}.
#'
#' @seealso \code{\link{daceBuilder}} \code{\link{daceLikelihood}} \code{\link{daceEvalFit}} 
#' @keywords internal
# @export
###################################################################################
daceGetFit <- function (theta, para){	
	if(any(theta <= 0)){ 
		stop('theta for dace model in spotPredictDace must be strictly positive')
	}
	
	fit = daceObjfunc(theta, para, "fit")$fit

	list(regr=para$regr, corr=para$corr, theta=theta, 
				beta=fit$beta, gamma=fit$gamma, sigma2=para$sY^2*fit$sigma2, 
				S=para$S, Ssc=para$Ssc, Y=para$y, Ysc=para$Ysc, C=fit$C, 
				Ft=fit$Ft, G=fit$G)#, detR=fit$detR);
}  
 
###################################################################################
#' Backslash operator.
#'
#' Reproduce what MATLAB's backslash operator can do, using qr() and qr.coef().
#'
#' @param X X matrix
#' @param Y Y vector
#'
#' @return Returns coefficients 
#' @keywords internal
# @export
###################################################################################
spotHelpBslash<-function(X,Y){
	X<-qr(X)
	qr.coef(X,Y)
}

###################################################################################
#' repmat
#'
#' Reproduce what MATLAB's repmat function can do, using kronecker function.
#'
#' @param a matrix
#' @param n dimension 1 (rows)
#' @param m dimension 2 (cols)
#'
#' @return Returns repeated matrix
#' @keywords internal
# @export
###################################################################################
repmat <- function(a,n,m) {kronecker(matrix(1,n,m),a)}

###################################################################################
#' DACE objective function
#'
#' Evaluate the fit of a certain set of model parameters (\code{theta}), and
#' get all relevant variables of the model.
#'
#' @param theta model parameters to be evaluated
#' @param para model option list, as created with \code{\link{dacePrepareFit}}
#' @param what a string: "all" both the likelihood (f) and the model list (fit) will be returned, "f" and "fit specify to return only those.
#'
#' @return A list of two elements (which are NA if \code{what} is specified accordingly)\cr
#' 	\code{f} likelihood \cr
#' 	\code{fit} also a list list of model variables, with the following elements: 
#' 		\item{\code{beta}}{ generalized least squares estimate}
#' 		\item{\code{gamma}}{ correlation factors}
#' 		\item{\code{sigma2}}{ maximum Likelihood estimate of the process variance}
#' 		\item{\code{C}}{ Cholesky factor of correlation matrix}
#' 		\item{\code{Ft}}{ Decorrelated regression matrix}
#' 		\item{\code{G}}{ From QR factorization: Ft = Q*t(G)\cr}
#'
#' @author The authors of the original DACE Matlab code are Hans Bruun Nielsen \email{hbn@@imm.dtu.dk}, Soren Nymand Lophaven and Jacob Sondergaard. \cr
#' Extension of the Matlab code by Tobias Wagner \email{wagner@@isf.de}. \cr 
#' Porting and adaptation to R and further extensions by Martin Zaefferer \email{martin.zaefferer@@fh-koeln.de}.
#'
#' @seealso \code{\link{daceBuilder}} \code{\link{daceLikelihood}} \code{\link{daceEvalFit}} 
#' @keywords internal
# @export
###################################################################################
daceObjfunc <- function(theta, para, what="all"){
	#% Initialize
	obj = Inf #penalty for bad numerical conditions    TODO: Inf crashes some optimizers, some not. Use very large number instead?
	fit=NA
	m = dim(para$F)[1]
	if(is.null(m))m=length(para$F)
	#% Set up  R
	r = para$corr(theta, para$D, "r")$r
	idx = which(r>0);
	o = 1:m;   
	mu = (10+m)*.Machine$double.eps
	#R = sparseMatrix(
	#		i=c(para$ij[idx,1], o), 
	#		j=c(para$ij[idx,2], o),
	#		x=c(r[idx], rep(1,m)+mu)  #TODO? matrix(1,m,1) statt rep
	#	)	
	R<-matrix(0,nrow=m,ncol=m)   #using sparse matrixes is far to slow in R.
	R[cbind(c(para$ij[idx,1], o),c(para$ij[idx,2], o))]=c(r[idx], rep(1,m)+mu)
	#	browser()
	#% Cholesky factorization. If it fails, return Inf for f, and NA for fit
	#if(!is.positive.definite(R)) return(list(f=obj,fit=fit));     #remark: is.positive.definite(R) does nothing else but a try(chol(R)) if method chol is used. therefore skipped.
	Cmat = try( chol(R), TRUE)
	if(class(Cmat) == "try-error"){		
		return(list(f=obj,fit=fit))
	}
	
	
	#% Get least squares solution
	Cmat = t(Cmat)
	Ft = spotHelpBslash(Cmat,para$F)
	resqr = qr(Ft)
	Q=qr.Q(resqr)
	G=qr.R(resqr)
	
	if  (rcond(G) < 1e-10){
		#% Check   F  
		if (kappa(para$F,exact=T) > 1e15 ){ #TODO?
			stop('F is too ill conditioned. Poor combination of regression model and design sites');
		}else{  #% Matrix  Ft  is too ill conditioned
			return(list(f=obj,fit=fit))	
		} 
	}
	Yt = spotHelpBslash(Cmat,para$y)   
	beta = spotHelpBslash(G,(t(Q)%*%Yt))
	rho = if(length(beta)==1){Yt - Ft*beta}else{Yt - Ft%*%beta}
	sigma2 = colSums(rho^2)/m
	detR = prod( (diag(Cmat)) ^ (2/m) )
	if(what=="f" || what=="all") 
		obj = sum(sigma2) %*% detR
	if(what=="fit" || what=="all")
		fit = list(sigma2=sigma2, beta=beta, gamma= t(rho)%*%solve(Cmat), C=Cmat, Ft=Ft, G=t(G))
	list(f=obj,fit=fit)
}

###################################################################################
#' Correlation: Noisy Gauss
#'
#' Noisy Gaussian correlation function using nuggets.
#' 
#            n
#    r_i = prod exp(-theta_j * d_ij^2) ,  i = 1,...,m
#           j=1
#'
#' @param theta parameters in the correlation function
#' @param d m*n matrix with differences between given data points
#' @param ret A string. If set to \code{"all"} or \code{"dr"}, the derivative of \code{r} (\code{dr}) will be returned, else \code{dr} is \code{NA}.
#'
#' @return returns a list with two elements:
#' 			\item{\code{r}}{correlation}
#' 			\item{\code{dr}}{m*n matrix with the Jacobian of \code{r} at \code{x}. It is
#'           assumed that \code{x} is given implicitly by \code{d[i,] = x - S[i,]},
#'           where \code{S[i,]} is the \code{i}'th design site.}
#'
#' @seealso \code{\link{spotPredictDace}} \code{\link{daceBuilder}}
#'
#' @author The authors of the original DACE Matlab code \url{http://www2.imm.dtu.dk/~hbni/dace/} 
#' are Hans Bruun Nielsen \email{hbn@@imm.dtu.dk}, Soren Nymand Lophaven and Jacob Sondergaard. \cr
#' Extension of the Matlab code by Tobias Wagner \email{wagner@@isf.de}. \cr 
#' Ported to R by Martin Zaefferer \email{martin.zaefferer@@fh-koeln.de}.
#'
#' @export
#' @keywords internal
###################################################################################
corrnoisygauss = function(theta, d , ret="all"){
	siz = dim(d);  #% number of differences and dimension of data
	m = siz[1]
	n = siz[2]
	lt = length(theta);
	if(lt != (n + 1)){
		stop(paste('Length of theta must be',n+1)) 
	}

	pow = 2
	tt = matrix(-theta[1:n],m,n,byrow=TRUE) 
	nugget = theta[n+1];
	td = abs(d)^pow * tt;
	etd=exp(td)
	r=1
	for(i in 1 : ncol(td)){
		r=r*etd[,i]
	}
	r= nugget * r
	if(ret=="all" || ret=="dr"){
		dr = nugget * pow * tt * sign(d) * (abs(d)^(matrix(1,m,n))) * matrix(r,m,n,byrow=FALSE)#repmat(r,1,n);
	}else{	dr= NA	}
	list(r=r,dr=dr)
}

###################################################################################
#' Correlation: Gauss
#'
#' Gaussian correlation function, no nugget.\cr
#' If \code{length(theta) = 1}, then the model is isotropic:\cr
#' all \code{theta_j = theta}.
#           n
#   r_i = prod exp(-theta_j * d_ij^2) ,  i = 1,...,m
#          j=1
#'
#' @param theta parameters in the correlation function
#' @param d m*n matrix with differences between given data points
#' @param ret A string. If set to \code{"all"} or \code{"dr"}, the derivative of \code{r} (\code{dr}) will be returned, else \code{dr} is \code{NA}.
#'
#' @return returns a list with two elements:
#' 			\item{\code{r}}{correlation}
#' 			\item{\code{dr}}{m*n matrix with the Jacobian of \code{r} at \code{x}. It is
#'           assumed that \code{x} is given implicitly by \code{d[i,] = x - S[i,]},
#'           where \code{S[i,]} is the \code{i}'th design site.}
#'
#' @seealso \code{\link{spotPredictDace}} \code{\link{daceBuilder}}
#'
#' @author The authors of the original DACE Matlab code \url{http://www2.imm.dtu.dk/~hbni/dace/} 
#' are Hans Bruun Nielsen \email{hbn@@imm.dtu.dk}, Soren Nymand Lophaven and Jacob Sondergaard. \cr
#' Ported to R by Martin Zaefferer \email{martin.zaefferer@@fh-koeln.de}.
#'
#' @export
#' @keywords internal
###################################################################################
corrgauss <- function(theta,d,ret="all"){
	siz = dim(d);  #% number of differences and dimension of data
	m = siz[1]
	n = siz[2]
	lt = length(theta);
	if (lt == 1){
		theta = rep(theta,n)
	}else if(lt != n){
		stop(paste("Length of theta must be ",n))
	}
	pow = 2
	tt = matrix(-theta[1:n],m,n,byrow=TRUE)  #TODO vllt t() hier und oben weglassen?
	td = d^pow * tt;
	r=exp(rowSums(td))

	if(ret=="all" || ret=="dr"){
		dr = matrix(-2*theta[1:n],m,n,byrow=TRUE) * d * matrix(r,m,n,byrow=FALSE)
	}else{	dr= NA	}
	list(r=r,dr=dr)
}

#TODO docu: equations?
###################################################################################
#' Correlation: Noisy Kriging
#'
#' Noisy Kriging correlation function using nuggets
#' 
#           n
#   r_i = prod exp(-theta_j * d_ij^theta_n+j)
#          j=1
#'
#' @param theta parameters in the correlation function
#' @param d m*n matrix with differences between given data points
#' @param ret A string. If set to \code{"all"} or \code{"dr"}, the derivative of \code{r} (\code{dr}) will be returned, else \code{dr} is \code{NA}.
#'
#' @return returns a list with two elements:
#' 			\item{\code{r}}{correlation}
#' 			\item{\code{dr}}{m*n matrix with the Jacobian of \code{r} at \code{x}. It is
#'           assumed that \code{x} is given implicitly by \code{d[i,] = x - S[i,]},
#'           where \code{S[i,]} is the \code{i}'th design site.}
#'
#' @seealso \code{\link{spotPredictDace}} \code{\link{daceBuilder}}
#'
#' @author The authors of the original DACE Matlab code \url{http://www2.imm.dtu.dk/~hbni/dace/} 
#' are Hans Bruun Nielsen \email{hbn@@imm.dtu.dk}, Soren Nymand Lophaven and Jacob Sondergaard. \cr
#' Extension of the Matlab code by Tobias Wagner \email{wagner@@isf.de}. \cr 
#' Ported to R by Martin Zaefferer \email{martin.zaefferer@@fh-koeln.de}.
#'
#' @export
#' @keywords internal
###################################################################################
corrnoisykriging = function(theta, d, ret="all"){
	siz = dim(d);  #% number of differences and dimension of data
	m = siz[1]
	n = siz[2]
	lt = length(theta);
	if(lt != (2 * n + 1)){
		stop(paste('Length of theta must be',2*n+1)) 
	}

	pow = matrix(theta[(n+1):(2*n)],m,n,byrow=TRUE);
	#pow = theta[(n+1):(2*n)]
	tt =  matrix(-theta[1:n],m,n,byrow=TRUE)
	nugget = theta[2*n+1];
	td = abs(d)^pow * tt;


	etd=exp(td)
	r=1
	for(i in 1 : ncol(td)){
		r=r*etd[,i]
	}


	r= nugget * r
	if(ret=="all" || ret=="dr"){
		dr = nugget * pow * tt * sign(d) * (abs(d)^(pow - matrix(1,m,n))) * matrix(r,m,n,byrow=FALSE)#repmat(r,1,n);
	}else{	dr= NA	}
	list(r=r,dr=dr)
}


#TODO docu: equations?
###################################################################################
#' Correlation:  Kriging
#'
#' Kriging correlation function, no nugget
#' 
#           n
#   r_i = prod exp(-theta_j * d_ij^theta_n+j)
#          j=1
#'
#' @param theta parameters in the correlation function
#' @param d m*n matrix with differences between given data points
#' @param ret A string. If set to \code{"all"} or \code{"dr"}, the derivative of \code{r} (\code{dr}) will be returned, else \code{dr} is \code{NA}.
#'
#' @return returns a list with two elements:
#' 			\item{\code{r}}{correlation}
#' 			\item{\code{dr}}{m*n matrix with the Jacobian of \code{r} at \code{x}. It is
#'           assumed that \code{x} is given implicitly by \code{d[i,] = x - S[i,]},
#'           where \code{S[i,]} is the \code{i}'th design site.}
#'
#' @seealso \code{\link{spotPredictDace}} \code{\link{daceBuilder}}
#'
#' @author The authors of the original DACE Matlab code \url{http://www2.imm.dtu.dk/~hbni/dace/} 
#' are Hans Bruun Nielsen \email{hbn@@imm.dtu.dk}, Soren Nymand Lophaven and Jacob Sondergaard. \cr
#' Extension of the Matlab code by Tobias Wagner \email{wagner@@isf.de}. \cr 
#' Ported to R by Martin Zaefferer \email{martin.zaefferer@@fh-koeln.de}.
#'
#' @export
#' @keywords internal
###################################################################################
corrkriging <- function(theta,d,ret="all"){
	siz = dim(d);  #% number of differences and dimension of data
	m = siz[1]
	n = siz[2]
	lt = length(theta);
	if(lt != (2 * n)){
		stop(paste('Length of theta must be',2*n)) 
	}
	pow = matrix(theta[(n+1):(2*n)],m,n,byrow=TRUE);
	tt = matrix(-theta[1:n],m,n,byrow=TRUE)  #TODO vllt t() hier und oben weglassen?
	td = abs(d)^pow * tt;
	r=exp(rowSums(td))

	if(ret=="all" || ret=="dr"){
		dr = matrix(-2*theta[1:n],m,n,byrow=TRUE) * d * matrix(r,m,n,byrow=FALSE)
	}else{	dr= NA	}
	list(r=r,dr=dr)
}

###################################################################################
#' Correlation: Cubic
#'
#' Cubic correlation function.\cr
#' If \code{length(theta) = 1}, then the model is isotropic:\cr
#' all \code{theta_j = theta}.
#           n
#   r_i = prod max(0, 1 - 3(theta_j*d_ij)^2 + 2(theta_j*d_ij)^3) ,  i = 1,...,m
#          j=1
#'
#' @param theta parameters in the correlation function
#' @param d m*n matrix with differences between given data points
#' @param ret A string. If set to \code{"all"} or \code{"dr"}, the derivative of \code{r} (\code{dr}) will be returned, else \code{dr} is \code{NA}.
#'
#' @return returns a list with two elements:
#' 			\item{\code{r}}{correlation}
#' 			\item{\code{dr}}{m*n matrix with the Jacobian of \code{r} at \code{x}. It is
#'           assumed that \code{x} is given implicitly by \code{d[i,] = x - S[i,]},
#'           where \code{S[i,]} is the \code{i}'th design site.}
#'
#' @seealso \code{\link{spotPredictDace}} \code{\link{daceBuilder}}
#'
#' @author The authors of the original DACE Matlab code \url{http://www2.imm.dtu.dk/~hbni/dace/} 
#' are Hans Bruun Nielsen \email{hbn@@imm.dtu.dk}, Soren Nymand Lophaven and Jacob Sondergaard. \cr
#' Ported to R by Martin Zaefferer \email{martin.zaefferer@@fh-koeln.de}.
#'
#' @export
#' @keywords internal
###################################################################################
corrcubic <- function(theta,d , ret="all"){
	siz = dim(d);  #% number of differences and dimension of data
	m = siz[1]
	n = siz[2]
	lt = length(theta);
	if (lt == 1){
		theta = rep(theta,n)
	}else if(lt != n){
		stop(paste('Length of theta must be 1 or ',n)) 
	}#else
	  #theta = theta(:).'; #force line vector, should be unnecessary in R
	#end
	#td = min(abs(d) * matrix(theta,m,lt,byrow=TRUE), 1);
	td = pmin(abs(d) * matrix(theta,m,lt,byrow=TRUE),1)
	ss = 1 - td^2 * (3 - 2*td);
	if(is.matrix(ss)){
		r = apply(ss,1,prod)
	}else{
		r = prod(ss)
	}
	if(ret=="all" || ret=="dr"){
		#dr = matrix(-2*theta[1:n],m,n,byrow=TRUE) * d * matrix(r,m,n,byrow=FALSE)		
		dr = matrix(0,m,n);
		for (j in 1:n){
			dd = 6*theta[j] * sign(d[,j]) * td[,j] * (td[,j] - 1);
			st<-ss[,c(1:(j-1), (j+1):n)]
			if(is.matrix(st)){
				dr[,j] = apply(st[,c(1:(j-1), (j+1):n)],1,prod)* dd;
			}else{
				dr[,j] = prod(st[,c(1:(j-1), (j+1):n)]) * dd;
			}
		}		
	}else{	dr= NA	}
	list(r=r,dr=dr)
}

###################################################################################
#' Correlation: Exp
#'
#' Exponential correlation function.\cr
#' If \code{length(theta) = 1}, then the model is isotropic:\cr
#' all \code{theta_j = theta}.
#           n
#   r_i = prod exp(-theta_j * |d_ij|)
#          j=1
#'
#' @param theta parameters in the correlation function
#' @param d m*n matrix with differences between given data points
#' @param ret A string. If set to \code{"all"} or \code{"dr"}, the derivative of \code{r} (\code{dr}) will be returned, else \code{dr} is \code{NA}.
#'
#' @return returns a list with two elements:
#' 			\item{\code{r}}{correlation}
#' 			\item{\code{dr}}{m*n matrix with the Jacobian of \code{r} at \code{x}. It is
#'           assumed that \code{x} is given implicitly by \code{d[i,] = x - S[i,]},
#'           where \code{S[i,]} is the \code{i}'th design site.}
#'
#' @seealso \code{\link{spotPredictDace}} \code{\link{daceBuilder}}
#'
#' @author The authors of the original DACE Matlab code \url{http://www2.imm.dtu.dk/~hbni/dace/} 
#' are Hans Bruun Nielsen \email{hbn@@imm.dtu.dk}, Soren Nymand Lophaven and Jacob Sondergaard. \cr
#' Ported to R by Martin Zaefferer \email{martin.zaefferer@@fh-koeln.de}.
#'
#' @export
#' @keywords internal
###################################################################################
correxp <- function(theta,d,ret="all"){
	siz = dim(d);  #% number of differences and dimension of data
	m = siz[1]
	n = siz[2]
	lt = length(theta);
	if (lt == 1){
		theta = rep(theta,n)
	}else if(lt != n){
		stop(paste('Length of theta must be 1 or ',n)) 
	}
	td = abs(d) * matrix(-theta,m,n,byrow=TRUE) 
	r=exp(rowSums(td))

	if(ret=="all" || ret=="dr"){
		dr = matrix(-theta[1:n],m,n,byrow=TRUE) * sign(d) * matrix(r,m,n,byrow=FALSE)
	}else{	dr= NA	}
	list(r=r,dr=dr)
}

###################################################################################
#' Correlation: Expg
#'
#' General exponential correlation function.\cr
#' If \code{n > 1}  and \code{length(theta) = 1}, then the model is isotropic:\cr
#' \code{theta_j = theta[1], j=1,...,n;  theta_[n+1] = theta[2]}.
#            n
#    r_i = prod exp(-theta_j * d_ij^theta_n+1)
#           j=1
#'
#' @param theta parameters in the correlation function
#' @param d m*n matrix with differences between given data points
#' @param ret A string. If set to \code{"all"} or \code{"dr"}, the derivative of \code{r} (\code{dr}) will be returned, else \code{dr} is \code{NA}.
#'
#' @return returns a list with two elements:
#' 			\item{\code{r}}{correlation}
#' 			\item{\code{dr}}{m*n matrix with the Jacobian of \code{r} at \code{x}. It is
#'           assumed that \code{x} is given implicitly by \code{d[i,] = x - S[i,]},
#'           where \code{S[i,]} is the \code{i}'th design site.}
#'
#' @seealso \code{\link{spotPredictDace}} \code{\link{daceBuilder}}
#'
#' @author The authors of the original DACE Matlab code \url{http://www2.imm.dtu.dk/~hbni/dace/} 
#' are Hans Bruun Nielsen \email{hbn@@imm.dtu.dk}, Soren Nymand Lophaven and Jacob Sondergaard. \cr
#' Ported to R by Martin Zaefferer \email{martin.zaefferer@@fh-koeln.de}.
#'
#' @export
#' @keywords internal
###################################################################################
correxpg <- function(theta,d,ret="all"){

	siz = dim(d);  #% number of differences and dimension of data
	m = siz[1]
	n = siz[2]
	lt = length(theta);
	if (n > 1 && lt == 2){
		theta = c(rep(theta[1],n),theta[2])
	}else if(lt != n+1){
		stop(paste('Length of theta must be 2 or ',n+1)) 
	}
	
	pow = theta[n+1];
	tt = matrix(-theta[1:n],m,n,byrow=TRUE)  #TODO vllt t() hier und oben weglassen?
	td = abs(d)^pow * tt;
	r=exp(rowSums(td))

	if(ret=="all" || ret=="dr"){
		dr = pow * tt * sign(d) * (abs(d)^(pow-1)) * matrix(r,m,n,byrow=FALSE)
	}else{	dr= NA	}
	list(r=r,dr=dr)
}


###################################################################################
#' Correlation: Lin
#'
#' Linear correlation function.\cr
#' If \code{length(theta) = 1}, then the model is isotropic:\cr
#' all \code{theta_j = theta}.
#            n
#    r_i = prod max(0, 1 - theta_j * d_ij) ,  i = 1,...,m
#           j=1
#'
#' @param theta parameters in the correlation function
#' @param d m*n matrix with differences between given data points
#' @param ret A string. If set to \code{"all"} or \code{"dr"}, the derivative of \code{r} (\code{dr}) will be returned, else \code{dr} is \code{NA}.
#'
#' @return returns a list with two elements:
#' 			\item{\code{r}}{correlation}
#' 			\item{\code{dr}}{m*n matrix with the Jacobian of \code{r} at \code{x}. It is
#'           assumed that \code{x} is given implicitly by \code{d[i,] = x - S[i,]},
#'           where \code{S[i,]} is the \code{i}'th design site.}
#'
#' @seealso \code{\link{spotPredictDace}} \code{\link{daceBuilder}}
#'
#' @author The authors of the original DACE Matlab code \url{http://www2.imm.dtu.dk/~hbni/dace/} 
#' are Hans Bruun Nielsen \email{hbn@@imm.dtu.dk}, Soren Nymand Lophaven and Jacob Sondergaard. \cr
#' Ported to R by Martin Zaefferer \email{martin.zaefferer@@fh-koeln.de}.
#'
#' @export
#' @keywords internal
###################################################################################
corrlin <- function(theta,d , ret="all"){
	siz = dim(d);  #% number of differences and dimension of data
	m = siz[1]
	n = siz[2]
	lt = length(theta);
	if (lt == 1){
		theta = rep(theta,n)
	}else if(lt != n){
		stop(paste('Length of theta must be 1 or ',n)) 
	}
	td = pmax(1- abs(d) * matrix(theta,m,lt,byrow=TRUE),0)
	if(is.matrix(td)){
		r = apply(td,1,prod)
	}else{
		r = prod(td)
	}
	if(ret=="all" || ret=="dr"){
		dr = matrix(0,m,n);
		for (j in 1:n){
			dd = (-theta[j] * sign(d[,j]))
			st<-td[,c(1:(j-1), (j+1):n)]
			if(is.matrix(st)){
				dr[,j] = apply(st[,c(1:(j-1), (j+1):n)],1,prod)* dd;
			}else{
				dr[,j] = prod(st[,c(1:(j-1), (j+1):n)]) * dd;
			}
		}		
	}else{	dr= NA	}
	list(r=r,dr=dr)
}


###################################################################################
#' Correlation: Spherical
#'
#' Spherical correlation function.\cr
#' If \code{length(theta) = 1}, then the model is isotropic:\cr
#' all \code{theta_j = theta}.
#            n
#    r_i = prod max(0, 1 - 1.5(theta_j*d_ij) + .5(theta_j*d_ij)^3) ,  i = 1,...,m
#           j=1
#'
#' @param theta parameters in the correlation function
#' @param d m*n matrix with differences between given data points
#' @param ret A string. If set to \code{"all"} or \code{"dr"}, the derivative of \code{r} (\code{dr}) will be returned, else \code{dr} is \code{NA}.
#'
#' @return returns a list with two elements:
#' 			\item{\code{r}}{correlation}
#' 			\item{\code{dr}}{m*n matrix with the Jacobian of \code{r} at \code{x}. It is
#'           assumed that \code{x} is given implicitly by \code{d[i,] = x - S[i,]},
#'           where \code{S[i,]} is the \code{i}'th design site.}
#'
#' @seealso \code{\link{spotPredictDace}} \code{\link{daceBuilder}}
#'
#' @author The authors of the original DACE Matlab code \url{http://www2.imm.dtu.dk/~hbni/dace/} 
#' are Hans Bruun Nielsen \email{hbn@@imm.dtu.dk}, Soren Nymand Lophaven and Jacob Sondergaard. \cr
#' Ported to R by Martin Zaefferer \email{martin.zaefferer@@fh-koeln.de}.
#'
#' @export
#' @keywords internal
###################################################################################
corrspherical <- function(theta,d , ret="all"){
	siz = dim(d);  #% number of differences and dimension of data
	m = siz[1]
	n = siz[2]
	lt = length(theta);
	if (lt == 1){
		theta = rep(theta,n)
	}else if(lt != n){
		stop(paste('Length of theta must be 1 or ',n)) 
	}#else
	  #theta = theta(:).'; #force line vector, should be unnecessary in R
	#end
	#td = min(abs(d) * matrix(theta,m,lt,byrow=TRUE), 1);
	td = pmin(abs(d) * matrix(theta,m,lt,byrow=TRUE),1)
	ss = 1 - td * (1.5 - 0.5*td^2);
	if(is.matrix(ss)){
		r = apply(ss,1,prod)
	}else{
		r = prod(ss)
	}
	if(ret=="all" || ret=="dr"){
		#dr = matrix(-2*theta[1:n],m,n,byrow=TRUE) * d * matrix(r,m,n,byrow=FALSE)		
		dr = matrix(0,m,n);
		for (j in 1:n){
			dd = 1.5*theta[j] * sign(d[,j]) * (td[,j]^2 - 1);
			st<-ss[,c(1:(j-1), (j+1):n)]
			if(is.matrix(st)){
				dr[,j] = apply(st[,c(1:(j-1), (j+1):n)],1,prod)* dd;
			}else{
				dr[,j] = prod(st[,c(1:(j-1), (j+1):n)]) * dd;
			}
		}		
	}else{	dr= NA	}
	list(r=r,dr=dr)
}

###################################################################################
#' Correlation: Spline
#'
#' Cubic spline correlation function.\cr
#' If \code{length(theta) = 1}, then the model is isotropic:\cr
#' all \code{theta_j = theta}.
#            n
#    r_i = prod S(theta_j*d_ij) ,  i = 1,...,m
#           j=1
# 
#  with
#            1 - 15x^2 + 30x^3   for   0 <= x <= 0.5
#    S(x) =  1.25(1 - x)^3       for  0.5 < x < 1
#            0                   for    x >= 1
#'
#' @param theta parameters in the correlation function
#' @param d m*n matrix with differences between given data points
#' @param ret A string. If set to \code{"all"} or \code{"dr"}, the derivative of \code{r} (\code{dr}) will be returned, else \code{dr} is \code{NA}.
#'
#' @return returns a list with two elements:
#' 			\item{\code{r}}{correlation}
#' 			\item{\code{dr}}{m*n matrix with the Jacobian of \code{r} at \code{x}. It is
#'           assumed that \code{x} is given implicitly by \code{d[i,] = x - S[i,]},
#'           where \code{S[i,]} is the \code{i}'th design site.}
#'
#' @seealso \code{\link{spotPredictDace}} \code{\link{daceBuilder}}
#'
#' @author The authors of the original DACE Matlab code \url{http://www2.imm.dtu.dk/~hbni/dace/} 
#' are Hans Bruun Nielsen \email{hbn@@imm.dtu.dk}, Soren Nymand Lophaven and Jacob Sondergaard. \cr
#' Ported to R by Martin Zaefferer \email{martin.zaefferer@@fh-koeln.de}.
#'
#' @export
#' @keywords internal
###################################################################################
corrspline <- function(theta,d , ret="all"){
	siz = dim(d);  #% number of differences and dimension of data
	m = siz[1]
	n = siz[2]
	lt = length(theta);
	if (lt == 1){
		theta = rep(theta,n)
	}else if(lt != n){
		stop(paste('Length of theta must be 1 or ',n)) 
	}
	mn = m*n;   
	ss = matrix(0,mn,1);
	xi = matrix(abs(d) * matrix(theta,m,lt,byrow=TRUE), mn,1);
	#% Contributions to first and second part of spline
	i1 = which(xi <= 0.2);
	i2 = which(0.2 < xi & xi < 1);
	#if(length(i1)<1){
		ss[i1] = 1 - xi[i1]^2 * (15  - 30*xi[i1]);
	#}
	#if(length(i2)<1){
		ss[i2] = 1.25 * (1 - xi[i2])^3;
	#}
	#% Values of correlation
	ss = matrix(ss,m,n);
	if(is.matrix(ss)){
		r = apply(ss,1,prod)
	}else{
		r = prod(ss)
	}
 

	if(ret=="all" || ret=="dr"){ #% get Jacobian
		u = matrix(sign(d) * matrix(theta,m,lt,byrow=TRUE), mn,1);
		dr = matrix(0,mn,1);
		#if  ~isempty(i1)
			dr[i1,] = u[i1] * ( (90*xi[i1] - 30) * xi[i1] );
		#end
		#if  ~isempty(i2)
			dr[i2,] = -3.75 * u[i2] * (1 - xi[i2])^2;
		#end
		ii = 1 : m;
		for (j in 1:n){
			sj = ss[,j];  
			ss[,j] = dr[ii];
			if(is.matrix(ss)){
				dr[ii,]  = apply(ss,1,prod)
			}else{
				dr[ii,]  = prod(ss)
			}
			ss[,j] = sj;   
			ii = ii + m;
		}
		dr = matrix(dr,m,n);	
	}else{	dr= NA	}
	list(r=r,dr=dr)
}

###################################################################################
#' Regression: Regpoly0
#'
#' Zero order polynomial regression function.
#'
#' @param S m*n matrix with \code{m} design points of dimension \code{n}
#' @param grad define if function returns gradient, default is \code{FALSE}
#'
#' @return returns a list with two elements:
#' 			\item{\code{f}}{vector of ones, with length \code{m}}
#' 			\item{\code{df}}{Jacobian at the first point (first row in S) }
#'
#' @seealso \code{\link{spotPredictDace}} \code{\link{daceBuilder}}
#'
#' @author The authors of the original DACE Matlab code \url{http://www2.imm.dtu.dk/~hbni/dace/} 
#' are Hans Bruun Nielsen \email{hbn@@imm.dtu.dk}, Soren Nymand Lophaven and Jacob Sondergaard. \cr
#' Ported to R by Martin Zaefferer \email{martin.zaefferer@@fh-koeln.de}.
#'
#' @export
#' @keywords internal
###################################################################################
regpoly0 = function(S,grad=FALSE){
	siz = dim(S);  
	m = siz[1]
	n = siz[2]
	f = rep(1,m);
	
	if(grad){
		df = rep(0,n);
	}else{
		df = NULL
	}	
	
	list(f=f,df=df)
}

##################################################################################
#' Regression: Regpoly1
#'
#' First order polynomial regression function.
#'
#' @param S m*n matrix with \code{m} design points of dimension \code{n}
#' @param grad define if function returns gradient, default is \code{FALSE}
#'
#' @return returns a list with two elements:
#' 			\item{\code{f}}{matrix of two columns: \cr 
#'			1. vector of ones with length \code{m} \cr
#'			2. \code{S}}
#' 			\item{\code{df}}{Jacobian at the first point (first row in S) }
#'
#' @seealso \code{\link{spotPredictDace}} \code{\link{daceBuilder}}
#'
#' @author The authors of the original DACE Matlab code \url{http://www2.imm.dtu.dk/~hbni/dace/} 
#' are Hans Bruun Nielsen \email{hbn@@imm.dtu.dk}, Soren Nymand Lophaven and Jacob Sondergaard. \cr
#' Ported to R by Martin Zaefferer \email{martin.zaefferer@@fh-koeln.de}.
#'
#' @export
#' @keywords internal
###################################################################################
regpoly1 <- function(S,grad=FALSE){
	siz = dim(S);  
	m = siz[1]
	n = siz[2]
	f = cbind(rep(1,m), S);
	
	
	if(grad){
		df = cbind(rep(0,n), diag(rep(1,n)));
	}else{
		df = NULL
	}	
	
	list(f=f,df=df)
}


##################################################################################
#' Regression: Regpoly2
#'
#' Second order polynomial regression function.
#'
#' @param S m*n matrix with \code{m} design points of dimension \code{n}
#' @param grad define if function returns gradient, default is \code{FALSE}
#'
#' @return returns a list with two elements:
#' 			\item{\code{f}}{matrix: \cr 
#'			( 1 S S[,1]*S S[,2]S[,2:n] ... S[,n]^2 )}
#' 			\item{\code{df}}{Jacobian at the first point (first row in S) }
#'
#' @seealso \code{\link{spotPredictDace}} \code{\link{daceBuilder}}
#'
#' @author The authors of the original DACE Matlab code \url{http://www2.imm.dtu.dk/~hbni/dace/} 
#' are Hans Bruun Nielsen \email{hbn@@imm.dtu.dk}, Soren Nymand Lophaven and Jacob Sondergaard. \cr
#' Ported to R by Martin Zaefferer \email{martin.zaefferer@@fh-koeln.de}.
#'
#' @export
#' @keywords internal
###################################################################################
regpoly2 <- function(S,grad=FALSE){
	siz = dim(S);  
	m = siz[1]
	n = siz[2]
	nn = (n+1)*(n+2)/2;  #% Number of columns in f  
	#% Compute  f	
	f = cbind(rep(1,m), S,  matrix(0,m,nn-n-1));
	j = n+1;   
	q = n;
	for(k in 1:n){
	  f[,j+(1:q)] = matrix(S[,k],m,q) * S[,k:n];
	  j = j+q;   
	  q = q-1;
	}	
	if(grad){
		df = cbind(rep(0,n), diag(rep(1,n)),  matrix(0,n,nn-n-1));
		j = n+1;   
		q = n; 
		for(k in 1:n){
			if(k<n){ 
				df[k,j+(1:q)] = cbind(2*S[1,k], S[1,(k+1):n]);
				for(i in 1:(n-k)){
					df[k+i,j+1+i] = S[1,k]; 
				}
			}else{ df[k,j+(1:q)] = cbind(2*S[1,k]);}

			j = j+q;   
			q = q-1;
		}			
	}else{
		df = NULL
	}	
	list(f=f,df=df)
}

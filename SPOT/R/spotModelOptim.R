###################################################################################
#' Optimize predicted meta model
#'  
#' Optimizes an existing fit of a model to get an optimal new design point. 
#' Executed after building the prediction model in the sequential SPOT step.
#' 
#' @param startPoint initial points for the optimization
#' @param spotConfig list of all options, needed to provide data for calling functions.
#' 
#' This function uses the parameter \code{spotConfig$seq.modelFit}. This is supposed to be a fit, that can be evaluated by the associated \code{spotConfig$seq.predictionModel.func} function.
#' The parameter \code{spotConfig$seq.predictionOpt.method} will be used to choose the optimization method to be used to find the minimum of the fitted model:\cr
#' "optim-L-BFGS-B" - BFGS quasi-Newton: \code{stats} Package\cr
#' "pso" - Particle Swarm Optimization: \code{pso} Package \cr
#' "cmaes" - Covariance Matrix Adaptation Evolution Strategy: \code{cmaes} Package\cr
#' "genoud" - Combines evolutionary search algorithms with derivative-based (Newton or quasi-Newton) methods: \code{rgenoud} Package\cr
#' "DEoptim" - Differential Evolution implementation: \code{DEoptim} Package\cr
#' "bobyqa" - Trust region method that forms quadratic models by interpolation: \code{minqa} Package\cr
#' "BBoptim" - Strategy using different Barzilai-Borwein step-lengths: \code{BB} Package\cr
#' "GenSA" - Generalized simulated annealing which for global minimization of a very complex non-linear objective function with a very large number of optima: \code{GenSA} Package\cr
#' "hjkb" - Bounded Hooke-Jeeves algorithm for derivative-free optimization: \code{dfoptim} Package\cr\cr
#' Additionally to the above methods, several methods from the package \code{nloptr} can be chosen. For instance:\cr
#' "NLOPT_LN_NELDERMEAD" - Nelder-Mead Simplex\cr
#' "NLOPT_LN_SBPLX" - Nelder-Mead Simplex on sequence of subspaces\cr
#' "NLOPT_GN_DIRECT"  - Direct Search \cr
#' "NLOPT_GN_DIRECT_L"  - Direct Search, locally biased\cr\cr
#' The complete list of suitable nlopt methods (non-gradient, bound constraints) is: \cr
#'	"NLOPT_GN_DIRECT","NLOPT_GN_DIRECT_L","NLOPT_GN_DIRECT_L_RAND",
#'	"NLOPT_GN_DIRECT_NOSCAL","NLOPT_GN_DIRECT_L_NOSCAL","NLOPT_GN_DIRECT_L_RAND_NOSCAL",
#'	"NLOPT_GN_ORIG_DIRECT","NLOPT_GN_ORIG_DIRECT_L","NLOPT_LN_PRAXIS",							
#'	"NLOPT_GN_CRS2_LM","NLOPT_LN_COBYLA","NLOPT_LN_NEWUOA_BOUND",
#'	"NLOPT_LN_NELDERMEAD","NLOPT_LN_SBPLX","NLOPT_LN_BOBYQA","NLOPT_GN_ISRES"\cr\cr
#' All of the above methods use bound constraints, which will be chosen with the limits specified in \code{spotConfig$alg.roi} (or the \code{.roi} file).
#' For references and details on the specific methods, please check the documentation of the packages that provide them.
#' 
#' Note that some methods may require additional parametrization. For this purpose, it would be recommended to use \code{spotPredictOptMulti}
#' as a template to write custom functions. \code{spotPredictOptMulti} itself is limited to calling mostly default settings of a large number of different optimizers available in \code{R}.
#' 
#' @return returns the list \code{spotConfig} with two new entries:\cr
#' 	\code{spotConfig$optDesign} are the parameters of the new minimal design point \cr
#'	\code{spotConfig$optDesignY} is the associated value of the objective function
#' @seealso  \code{\link{spotModelParetoOptim}} solves the same task for multi-criteria optimization (i.e. more than just one surrogate model)\cr
#' @export
###################################################################################
spotModelOptim <- function(startPoint,spotConfig){
	startPoint=as.numeric(startPoint[which.min(spotConfig$seq.largeDesignY[[1]]),])
	spotWriteLines(spotConfig$io.verbosity,3,"spotModelOptim started")
	if(is.null(spotConfig$seq.predictionOpt.method)) spotConfig$seq.predictionOpt.method = "optim-L-BFGS-B"
	if(is.null(spotConfig$seq.predictionOpt.budget)) spotConfig$seq.predictionOpt.budget = 200
	if(is.null(spotConfig$seq.predictionOpt.psize)) spotConfig$seq.predictionOpt.psize = 20
	if(is.null(spotConfig$seq.predictionOpt.restarts)) spotConfig$seq.predictionOpt.restarts = FALSE
	fit<-spotConfig$seq.modelFit
	tempVerbose<-spotConfig$io.verbosity
	#building the local objective function, differences depend on type of model.
	spotGetFitness <- function(x){
		spotConfig$io.verbosity=0
		eval(call(spotConfig$seq.predictionModel.func
								, NULL 
								, NULL
								, as.data.frame(x)
								, spotConfig
								, fit #external fit is used, model is only evaluated not build, therefore the NULLS are no problem
								))$seq.largeDesignY[[1]]
	}
	lowROI<-as.numeric(spotConfig$alg.aroi[[1]])
	upROI<-as.numeric(spotConfig$alg.aroi[[2]])

	res <- spotOptimizationInterface(par=startPoint,fn=spotGetFitness,gr=NULL,lower=lowROI,upper=upROI,method=spotConfig$seq.predictionOpt.method,
						control=list(
							popsize = spotConfig$seq.predictionOpt.psize
							,fevals = spotConfig$seq.predictionOpt.budget
							,ineq_constr = spotConfig$alg.constraint.ineq
							,vectorize=TRUE
							,verbosity=spotConfig$io.verbosity
							,restarts=spotConfig$seq.predictionOpt.restarts))
							

	newDesignPrediction <- res$value
	newDesign <- res$par

	spotPrint(spotConfig$io.verbosity,3,newDesign)
	spotWriteLines(spotConfig$io.verbosity,3,"spotModelOptim finished")
	
	spotConfig$optDesign<-newDesign
	spotConfig$optDesignY<-newDesignPrediction
	spotConfig$io.verbosity<-tempVerbose
	spotConfig
}

###################################################################################
#' DEPRECATED
#'  
#' see \code{\link{spotModelOptim}}
#'
#' @param startPoint see \code{\link{spotModelOptim}}
#' @param spotConfig see \code{\link{spotModelOptim}}
#' 
#' @export
#' @keywords internal
###################################################################################
spotPredictOptMulti <- function(startPoint,spotConfig){
	warning("spotPredictOptMulti is a deprecated function, please use spotModelOptim")
	spotModelOptim(startPoint,spotConfig)
}

###################################################################################
#' Steepest Descent on RSM (linear model)
#'  
#' Optimizes an existing fit of a linear model created by the rsm function. Uses steepest descent method and adaptation of ROI alternatingly.
#' 
#' @param startPoint not used here
#' @param spotConfig list of all options, needed to provide data for calling functions.
#' 
#' @return returns the list \code{spotConfig} with two new entries:\cr
#' 	\code{spotConfig$optDesign} are the parameters of the new minimal design point \cr
#'	\code{spotConfig$optDesignY} is the associated value of the objective function
#' @seealso  \code{\link{spotModelParetoOptim}}, \code{\link{spotModelOptim}}
#' @export
###################################################################################
spotModelDescentLm <- function(startPoint,spotConfig){
	#startPoint=as.numeric(startPoint[which.min(spotConfig$seq.largeDesignY[[1]]),])
	spotWriteLines(spotConfig$io.verbosity,3,"spotModelDescentLM started")
	pNames <- row.names(spotConfig$alg.roi)
	nParam <- length(pNames)
	fit<-spotConfig$seq.modelFit
	tempVerbose<-spotConfig$io.verbosity
	if(is.null(spotConfig$seq.predictionOpt.method)) spotConfig$seq.predictionOpt.method = "optim-L-BFGS-B"
	if(is.null(spotConfig$seq.predictionOpt.budget)) spotConfig$seq.predictionOpt.budget = 200
	if(is.null(spotConfig$seq.predictionOpt.psize)) spotConfig$seq.predictionOpt.psize = 20
	if(is.null(spotConfig$seq.predictionOpt.restarts)) spotConfig$seq.predictionOpt.restarts = FALSE
	#building the local objective function, differences depend on type of model.
	spotGetFitness <- function(x){			
		spotConfig$io.verbosity=0;
		eval(call(spotConfig$seq.predictionModel.func
								, NULL 
								, NULL
								, as.data.frame(x)
								, spotConfig
								, fit 
								))$seq.largeDesignY[[1]]
	}
	lowROI<-as.numeric(spotConfig$alg.aroi[[1]])
	upROI<-as.numeric(spotConfig$alg.aroi[[2]])
	
	#default parameters of steepest descent approach
	if(is.null(spotConfig$seq.useGradient)) spotConfig$seq.useGradient <- TRUE
	if(is.null(spotConfig$seq.useCanonicalPath)) spotConfig$seq.useCanonicalPath <- TRUE
	
	if ((( max(spotConfig$alg.currentResult$STEP) %% 2) == 0) | !spotConfig$seq.useAdaptiveRoi){ #every second step
		if ( (spotConfig$seq.useGradient == TRUE)){ #analytical gradient descent
			if (spotConfig$seq.useCanonicalPath == TRUE){ # start at saddle point in both directions (canonical path analysis)
				steepestDesc <-  as.data.frame(canonical.path(fit, descent=TRUE, dist = seq(-0.2,0.2, by = 0.1))[,2:eval(nParam+1)])
			}
			else{ # start at origin in one direction
				steepestDesc <-  as.data.frame(steepest(fit, descent=TRUE, dist = seq(0.1,0.5, by = 0.1))[,2:eval(nParam+1)])
			}
			## ensure feasibility			
			steepestDesc<- steepestDesc[apply(steepestDesc, 1, function(x) all(x > -1 & x < 1)), ]
			## use model optimization if steepest descent is empty
			if (nrow(steepestDesc)==0){
				#tmp <- spotOptimizationInterface(par=startPoint,fn=spotGetFitness,gr=NULL,lower=lowROI,upper=upROI,method=spotConfig$seq.predictionOpt.method,
				tmp <- spotOptimizationInterface(par=lowROI + (upROI - lowROI) / 2,fn=spotGetFitness,gr=NULL,lower=lowROI,upper=upROI,method=spotConfig$seq.predictionOpt.method,
						control=list(
							popsize = spotConfig$seq.predictionOpt.psize
							,fevals = spotConfig$seq.predictionOpt.budget
							,vectorize=TRUE
							,verbosity=spotConfig$io.verbosity
							,restarts=spotConfig$seq.predictionOpt.restarts))
				M <- matrix(tmp$par,ncol=length(lowROI))
				yM <- tmp$value			
			}else{
				steepestDesc <- code2val(steepestDesc, codings = fit$coding)			
				M<-steepestDesc;
				yM<-spotGetFitness(steepestDesc);
			}
		}
		else{ # numerical optimization of model
			#tmp <- spotOptimizationInterface(par=startPoint,fn=spotGetFitness,gr=NULL,lower=lowROI,upper=upROI,method=spotConfig$seq.predictionOpt.method,
			tmp <- spotOptimizationInterface(par=lowROI + (upROI - lowROI) / 2,fn=spotGetFitness,gr=NULL,lower=lowROI,upper=upROI,method=spotConfig$seq.predictionOpt.method,
						control=list(
							popsize = spotConfig$seq.predictionOpt.psize
							,fevals = spotConfig$seq.predictionOpt.budget
							,vectorize=TRUE
							,verbosity=spotConfig$io.verbosity
							,restarts=spotConfig$seq.predictionOpt.restarts))
			M <- matrix(tmp$par,ncol=length(lowROI))
			yM <- tmp$value
		}
	}
	else{ #recalibration phase
		## now we have the evaluated path of the steepeest descent
		## take the best point as the new center point		
		nbest <- nrow(spotConfig$alg.currentBest)
		xB <- spotConfig$alg.currentBest[nbest,c(pNames)]
		xB <- val2code(xB, fit$coding)
		ds <- min( rep(1,nParam) - abs(xB))	
		#browser()	
		print(ds)
		if (ds > 0.1){
			spotPrint(spotConfig$io.verbosity,2,"ds:")
			spotPrint(spotConfig$io.verbosity,2,ds)
			lower <- xB - ds
			upper <- xB + ds      
			A <- rbind(lower, upper)	
			rownames(A)=NULL
			colnames(A)=names(xB)
			A <- t(code2val(A, fit$coding))
			A <- data.frame(pNames,A,spotConfig$alg.roi$type)
			colnames(A) <- c("name","lower", "upper" , "type")
			if(spotConfig$spot.fileMode){ #write aroi file
				spotWriteAroi(A,spotConfig$io.verbosity,spotConfig$io.columnSep,spotConfig$io.aroiFileName)	
			}
			spotConfig$alg.aroi<-spotROI(A$lower,A$upper,type=A$type,varnames=A$name);
			spotWriteLines(spotConfig$io.verbosity,2,"AROI modified. Execution with continued in the adapted ROI.");
			## generate a new design 
			spotConfig$seq.useAdaptiveRoi <- TRUE
			M <- spotCreateDesignFrF2(spotConfig) #TODO, problem with FrF2 function, buggy.
			yM <- spotGetFitness(M)
		}
		else{
			spotConfig$seq.useAdaptiveRoi <- TRUE
			M <- spotCreateDesignLhs(spotConfig, noDesPoints=spotConfig$seq.design.new.size)
			yM <- spotGetFitness(M)
		}	
	}	
	#########################	
	spotPrint(spotConfig$io.verbosity,3,M)
	spotWriteLines(spotConfig$io.verbosity,3,"spotModelDescentLM finished")
	
	spotConfig$optDesign<-M
	spotConfig$optDesignY<-as.matrix(yM)
	spotConfig$io.verbosity<-tempVerbose
	spotConfig
}
#   Copyright (c) 2014-2015 by Martin Zaefferer, Cologne University of Applied Sciences

###################################################################################
#' Combinatorial Efficient Global Optimization
#' 
#' Model-based optimization for combinatorial or mixed problems. Based on measures of distance or dissimilarity.
#' 
#' @param x Optional initial design as a list. If NULL (default), \code{creationFunction} (in \code{control} list) is used to create initial design. 
#' If \code{x} has less individuals than specified by \code{control$evalInit}, \code{creationFunction} will fill up the design.
#' @param fun target function to be minimized
#' @param control (list), with the options of optimization and model building approaches employed\cr
#' \code{evalInit} Number of initial evaluations (i.e., size of the initial design), integer, default is \code{2}\cr
#' \code{vectorized} Boolean. Defines whether target function is vectorized (takes a list of solutions as argument) or not (takes single solution as argument). Default: FALSE\cr
#' \code{verbosity} Level of text output during run. Defaults to 0, no output.\cr
#' \code{plotting} Plot optimization progress during run (TRUE) or not (FALSE). Default is FALSE.\cr
#' \code{targetY} optimal value to be found, stopping criterion, default is \code{-Inf}\cr
#' \code{evalBudget} maximum number of target function evaluations, default is \code{100}\cr
#' \code{creationRetries} When a model does not predict an actually improving solution, a random exploration step is performed. \code{creationRetries} solutions are created randomly. 
#' 		For each, distance to all known solutions is calculated. The minimum distance is recorded for each random solution. 
#' 		The random solution with maximal minimum distance is chosen doe be evaluated in the next iteration.\cr
#' \code{model} Model to be used as a surrogate of the target function. Default is "K" (Kriging). Also
#'		available are: "LM" (linear, distance-based model), "RBFN" Radial Basis Function Network.\cr
#' \code{modelSettings} List of settings for model building, passed on as the control argument to the model training functions \code{\link{modelKriging}}, \code{\link{modelLinear}}, \code{\link{modelRBFN}}.\cr
#' \code{infill} This parameter specifies a function to be used for the infill criterion (e.g., the default is expected improvement \code{infillExpectedImprovement}).
#' To use no specific infill criterion this has to be set to \code{NA}. Infill criteria are only used with models that may provide some error estimate with predictions.\cr
#' \code{optimizer} Optimizer that finds the minimum of the surrogate model. Default is \code{"EA"} an Evolutionary Algorithm. No alternatives implemented yet.\cr
#' \code{optimizerSettings} List of settings for the method to optimize the model. \cr
#' \code{creationFunction} Function to create individuals/solutions in search space. Default is a function that creates random permutations of length 6\cr
#' \code{distanceFunction} distanceFunction a suitable distance function of type f(x1,x2), returning a scalar distance value, preferably between 0 and 1.
#'      Maximum distances larger 1 are not a problem, but may yield scaling bias when different measures are compared.
#' 		Should be non-negative and symmetric. With the setting \code{control$model="K"} this can also be a list of different fitness functions.
#'    Default is Hamming distance for permutations: distancePermutationHamming.
#'
#' @return a list:\cr 	
#' \code{xbest} best solution found\cr
#' \code{ybest} fitness of the best solution\cr
#' \code{x} history of all evaluated solutions\cr
#' \code{y} corresponding target function values f(x)\cr
#' \code{fit} model-fit created in the last iteration\cr
#' \code{fpred} prediction function created in the last iteration\cr
#' \code{count} number of performed target function evaluations\cr 
#' \code{message} message string, giving information on termination reason
#' \code{convergence} error/status code: \code{-1} for termination due 
#' to failed model building, \code{0} for termination due to depleted budget, 
#' \code{1} if attained objective value is equal to or below target (\code{control$targetY})
#' 
#' @examples
#' seed <- 0
#' glgseed=1
#' #distance
#' dF <- distancePermutationHamming
#' #mutation
#' mF <- mutationPermutationSwap
#' #recombination
#' rF <-  recombinationPermutationCycleCrossover 
#' #creation
#' cF <- function()sample(5)
#' #objective function
#' lF <- landscapeGeneratorUNI(1:5,dF)
#' #start optimization
#' set.seed(seed)
#' res1 <- optimCEGO(,lF,list(
#'				creationFunction=cF,
#'				distanceFunction=dF,
#'				optimizerSettings=list(budget=100,popsize=10,
#'				mutationFunction=mF,recombinationFunction=rF),
#'		evalInit=5,budget=15,targetY=0,verbosity=1,model=modelKriging,
#'		vectorized=TRUE)) ##target function is "vectorized", expects list as input
#' set.seed(seed)
#' res2 <- optimCEGO(,lF,list(
#'				creationFunction=cF,
#'				distanceFunction=dF,
#'				optimizerSettings=list(budget=100,popsize=10,
#'				mutationFunction=mF,recombinationFunction=rF),
#'				evalInit=5,budget=15,targetY=0,verbosity=1,model=modelRBFN,
#'		vectorized=TRUE)) ##target function is "vectorized", expects list as input
#' res1$xbest 
#' res2$xbest 
#'
#' @seealso \code{\link{modelKriging}}, \code{\link{modelLinear}}, \code{\link{modelRBFN}}, \code{\link{buildModel}}, \code{\link{optimEA}} 
#' 
#' @references Zaefferer, Martin; Stork, Joerg; Friese, Martina; Fischbach, Andreas; Naujoks, Boris; Bartz-Beielstein, Thomas. (2014). Efficient global optimization for combinatorial problems. In Proceedings of the 2014 conference on Genetic and evolutionary computation (GECCO '14). ACM, New York, NY, USA, 871-878. DOI=10.1145/2576768.2598282 http://doi.acm.org/10.1145/2576768.2598282 
#' @references Zaefferer, Martin; Stork, Joerg; Bartz-Beielstein, Thomas. (2014). Distance Measures for Permutations in Combinatorial Efficient Global Optimization. In Parallel Problem Solving from Nature - PPSN XIII (p. 373-383). Springer International Publishing.
#'
#' @export
###################################################################################
optimCEGO <- function(x=NULL,fun,control=list()){
	## default settings
	con<-list(evalInit = 2
			, vectorized = FALSE
			, verbosity = 0
			, plotting = FALSE
			, targetY = -Inf
      , budget = 100
			, creationRetries = 100
			, distanceFunction = distancePermutationHamming
			, creationFunction = solutionFunctionGeneratorPermutation(6)
			, infill= infillExpectedImprovement			
      , model = modelKriging
			, modelSettings= list()
      , optimizer = optimEA
			, optimizerSettings = list()
			, initialDesign = designRandom #no others available yet.
			, initialDesignSettings = list())
	con[names(control)] <- control
	control<-con
	rm(con)
	count <- control$evalInit
  vectorized <- control$vectorized
  verbosity <- control$verbosity
  plotting <- control$plotting
	creationFunction <- control$creationFunction
	distanceFunction <- control$distanceFunction
	
	## if target function not vectorized: vectorize with lapply
	fun #lazy load
	if(!vectorized) 
		fn <- function(x)unlist(lapply(x,fun))
	else 
		fn <- fun
	
	## Create main object of this function, which will also be the return value	
	res <- list(xbest=NA, ybest=NA, x=NA,y=NA,fit=NA,fpred=NA,distances=NA,count=count)
	
  ## Termination information:
  msg <- "Termination message:"
  res$convergence <- 0 #error/status code related to msg
  
	## Create initial design of experiment
	res$x <- control$initialDesign(x,creationFunction,count,control$initialDesignSettings)
	
	## Calculate distances between samples						TODO: what if distance function has parameters. see also kernel parameters. has to be calculated on the fly. switch off here?
	res$distances <- distanceMatrixWrapper(res$x,distanceFunction)
		
	## Evaluate initial population	
	res$y <- fn(res$x)

	## determine best
	indbest <- which.min(res$y)
	res$ybest <- res$y[[indbest]]
	res$xbest <- res$x[[indbest]]	
	
	## build initial model
	model <- buildModel(res,distanceFunction,control)

	## check whether EI infill is used
	useEI <- is.function(control$infill)
	
	## main loop
	while((res$count < control$budget) & (res$ybest > control$targetY)){
		## Optimize the surrogate:
    if(!identical(model,NA)){ 
      optimres <- optimizeModel(res,creationFunction,model,control)				
      ## Handle duplicated new candidate solution	
      duplicate <- optimres$xbest %in% res$x
			## Check whether next candidate solution is better than the best observed so far (based on model prediction)
      improved <- optimres$ybest < optimres$fpredbestKnownY
    }else{# model building failed, force termination
      msg <- paste(msg,"Model building failed, optimization stopped prematurely.")		
      warning("Model building failed in optimCEGO, optimization stopped prematurely.")
      res$convergence <- -1
      break;
    }
    
  	## Update evaluation counter
		res$count <- res$count+1
    
		if(!duplicate && ((improved || useEI))){ #exploitation step		#for a new individual to be accepted, it has to have a better predicted value than the prediction for the best known individual. One could also use the actual value of the best known individual, but this would deteriorate results in case of an offset.
			res$x[[res$count]] <- optimres$xbest		#NOTE; exploitation step and exploration is both the same when EI is used., thus the "||"
		}else{ #exploration step: no promising individual found, or promising individual is a duplicate -> create a new one randomly
			xc <- list() #candidates
			distlist <- NULL
			for (i in 1:control$creationRetries){ #create candidates, to get one new individual that is maximally different from the existing.
				xc[[i]] <- creationFunction()
				distx <- NULL
				if(length(distanceFunction)>1)
					dfn <- distanceFunction[[1]]
				else 
					dfn <- distanceFunction
				for(j in 1:length(res$x)){#calculate distance to each sampled location				
					distx <- c(distx,dfn(xc[[i]],res$x[[j]]))	#todo apply?
				}
				distlist <- c(distlist,min(distx)) #maximize the minimum distance
			}
			res$x[[res$count]] <- xc[[which.max(distlist)]] #this maximizes distance, but may still create a duplicate if max(min(dist))==0, e.g. if all randomly created individuals are duplicates of known solutions
		}
		res$x <- removeDuplicates(res$x, creationFunction)
		
		## evaluate with real target function
		res$y <-  c(res$y,fn(res$x[res$count]))
			
		## Logging
		indbest <- which.min(res$y)
		res$ybest <- res$y[[indbest]]
		res$xbest <- res$x[[indbest]]
		
		## Plotting and Text output
    if(verbosity > 0)
      print(paste("Evaluations:",res$count,"    Quality:",res$ybest))
		if(plotting){
			plot(res$y,type="l",xlab="number of evaluations", ylab="y")
			abline(res$ybest,0,lty=2)
		}

		## Update the distance matrix				#TODO what if distance parameters?
		res$distances <- distanceMatrixUpdate(res$distances,res$x,distanceFunction)
	
		## Update surrogate model and prediction function:
		model <- buildModel(res,distanceFunction,control)
	}
  
  #stopping criteria information for user:
	if(min(res$ybest,na.rm=TRUE) <= control$targetY) {
		msg <- paste(msg,"Successfully achieved target fitness.")		
    res$convergence <- 1
	}else if(res$count >= control$budget){ #budget exceeded
    msg <- paste(msg,"Target function evaluation budget depleted.")		
  }
  res$message <- msg
  
	res #return
}


###################################################################################
#' Model building
#' 
#' Model building support function for optimCEGO. Should not be called directly.
#' 
#' @param x list of samples in input space
#' @param y  matrix, column vector of observations for each sample
#' @param distanceFunction a suitable distance function of type f(x1,x2), returning a scalar distance value, preferably between 0 and 1.
#'      Maximum distances larger 1 are no problem, but may yield scaling bias when different measures are compared.
#' 		Should be non-negative and symmetric.  In case Kriging is chosen, it can also be a list of several distance functions. In this case, MLE is used 
#'		to determine the most suited distance measure (see the last reference).
#' @param control list of options \cr
#' \code{model} Model to be used as a surrogate of the target function. Default is "K" (Kriging). Also
#'		available are: "LM" (linear, distance-based model), "RBFN" Radial Basis Function Network.\cr
#' \code{modelSettings} List of settings for model building, passed on as the control argument to the model training functions \code{\link{modelKriging}}, \code{\link{modelLinear}}, \code{\link{modelRBFN}}.
#' \code{infill} This parameter specifies a function to be used for the infill criterion (e.g., the default is expected improvement \code{infillExpectedImprovement}).
#' To use no specific infill criterion this has to be set to \code{NA}. Infill criteria are only used with models that may provide some error estimate with predictions.\cr
#'
#' @return a list:\cr 	
#' \code{fit} model-fit \cr
#' \code{fpred} prediction function
#' 
#' @seealso \code{\link{optimCEGO}} 
#'
#' @keywords internal
#' @export
###################################################################################
buildModel <- function(res,distanceFunction,control){ 
	y <- res$y
	x <- res$x
	control$modelSettings$distances <- res$distances
	
	y #against lazy evaluation
	if(identical(control$model,"RBFN"))
		control$model <- modelRBFN
	else if(identical(control$model,"LM"))
		control$model <- modelLinear
	else if(identical(control$model,"K"))
		control$model <- modelKriging
	
	fit<-try(control$model(x,y,distanceFunction,control$modelSettings),TRUE)
	if(class(fit) == "try-error"){
    #warning("Model building in optimCEGO failed.") #same warning given in optimCEGO function
    return(NA)
  }
	
	fit$predAll <- is.function(control$infill)
	fit
	
	if(is.function(control$infill)){
		fpred <- function(x){	
				res=predict(fit,x)
				control$infill(res$y,res$s,min(y))
			} 
	}else{
		fpred <- function(x){predict(fit,x)$y}
	}

	list(fit=fit,fpred=fpred)
}


###################################################################################
#' Optimize Surrogate Model
#'
#' Interface to the optimization of the surrogate model
#' 
#' @param res TODO
#' @param creationFunction TODO
#' @param model TODO
#'
#' @return TODO
#' 
#' @seealso \code{\link{optimCEGO}} 
#'
#' @keywords internal
#' @export
###################################################################################
optimizeModel <- function(res,creationFunction,model,control){ 
	if(identical(control$optimizer,"EA")){
		control$optimizer=optimEA
	}
	if(is.null(control$optimizerSettings$creationFunction))
		control$optimizerSettings$creationFunction <- creationFunction
	optimres <- control$optimizer(NULL,model$fpred,control$optimizerSettings)
	optimres$fpredbestKnownY <- model$fpred(res$xbest) 
	optimres
}
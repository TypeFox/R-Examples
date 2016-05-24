### Feb, 1 2011: spot.noise added
### Jan, 29th 2010: alg.aroi added

###################################################################################################
#' spotOptim: optim-like spot interface
#'
#' Besides \code{\link{spot}} this is one of the main interfaces for using the SPOT package. It is build
#' like the \code{\link{optim}} interface.
#' 
#' It is of important to note that spot by default expects to optimize noisy functions. That means, the default settings of spot,
#' which are also used in spotOptim, include repeats of the initial and sequentially created design points. Also, as a default OCBA
#' is used to spread the design points for optimal usage of the function evaluation budget. OCBA will not work when there is no variance in the data.
#' So if the user wants to optimize non-noisy functions, the following settings should be used:\cr
#' \code{control$spot.ocba <- FALSE}\cr
#' \code{control$seq.design.maxRepeats <- 1}\cr
#' \code{control$init.design.repeats <- 1}\cr
#'
#' A call to a noisy function could look like this: \cr
#' \code{objFunction<-function(x){y=(x[1]+2)^2*(x[2]-4)^2+runif(1)}} \cr
#' \code{spotOptim(par=c(1,1),fn<-objFunction,lower=c(-10,-10),upper=c(10,10),method="spotPredictRandomForest",control=list(maxit=50))}
#'
#' A call to a non-noisy function could look like this: \cr
#' \code{objFunction<-function(x){y=(x[1]+2)^2*(x[2]-4)^2}} \cr
#' \code{spotOptim(par=c(1,1),fn<-objFunction,lower=c(-10,-10),upper=c(10,10),method="spotPredictRandomForest",control=list(maxit=50,spot.ocba=FALSE,seq.design.maxRepeats=1,init.design.repeats=1))}
#' 
#' @param par is a point in search interval (defines dimension)
#' @param fn is the target function (it can also be a string with the name of a spot interface function, like "\code{\link{spotFuncStartBranin}}")
#' @param gr gradient function, not implemented yet
#' @param ... additional parameters to be passed on to \code{fn}
#' @param lower is a vector that defines the lower boundary of search space
#' @param upper is a vector that defines the upper boundary of search space
#' @param method is a string that describes which method is to be used.
#' @param control is a list of additional settings. \code{maxit} is the number of function evaluations,
#' all the other settings will simply be passed to SPOT (see \code{\link{spotGetOptions}} for details) 
#'
#' @return This function returns a list with:\cr
#'	\code{par} parameters of the found solution\cr
#'	\code{value} target function value of the found solution\cr
#	\code{convergence} indicates successful completion when \code{0}\cr
#	\code{message}\cr
#	\code{hessian}\cr
#
#' @seealso \code{\link{spot}} \code{\link{spotOptimInterface}} \code{\link{spotOptimizationInterface}} \code{\link{spotOptimizationInterfaceMco}}
#'
#' @export
###################################################################################################
spotOptim <- function(par=NULL
                      , fn
                      , gr=NULL
					  , ...
                      , lower=-Inf
                      , upper=Inf
                      , method
                      ,control=list()){					  
	if (length(par)==0) stop("dimension of par is null")
	con<-list(maxit=100 #CON: Internal List with defaults for control
             , seed=123
             , spot.fileMode = FALSE);
	con[(namc <- names(control))] <- control;
	control<-con;
	
	control$alg.func<-fn; #will  be checked for string or function in spotGetOptions, here it is only passed
	control$seq.predictionModel.func<- method;
	if(is.null(control$auto.loop.nevals)){
		control$auto.loop.nevals<-control$maxit;
	}
	if(is.null(control$spot.seed)){
		control$spot.seed<- control$seed; 
	}

	# transform parameters for SPOT			
	### lower & upper ... region of interest
	control$alg.roi<-spotROI(lower,upper,dimROI=length(par))
	control$alg.aroi<-control$alg.roi;

	##### call SPOT	###############################			
	spotResult<-spot(spotConfig=control,...)
	
	##### prepare Results	###############################		
	result <- list()
	pNames <- row.names(spotResult$alg.roi);
	result$par <- as.numeric(spotResult$alg.currentBest[nrow(spotResult$alg.currentBest),pNames])
	result$value <- spotResult$alg.currentBest[nrow(spotResult$alg.currentBest),1]
	#result$counts <-
	result$convergence <- 0 # 0 indicates successful completion
	result$message <- NULL
	result$hessian <- NULL	
	result
}

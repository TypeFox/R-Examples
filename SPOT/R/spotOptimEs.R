
###################################################################################################
#' spotOptimEs: optim-like ES interface
#'
#' This is an interface to the Evolution Strategy used as a target algorithm by some SPOT demos. It is build
#' like the \code{\link{optim}} interface.
#' 
#' @param par is a point in search interval (defines dimension)
#' @param fn is the target function
#' @param gr gradient function, not used by this function
#' @param ... additional parameters to be passed on to \code{fn}
#' @param lower is a vector that defines the lower boundary of search space
#' @param upper is a vector that defines the upper boundary of search space
#' @param method this parameter is not used in the current version.
#' @param control is a list of additional settings. 
#'
#' The \code{control} list can contain the following settings:
#' \describe{
#' 			\item{maxit}{number of iterations, stopping criterion, default is \code{100}}
#' 			\item{mue}{number of parents, default is \code{10}}
#' 			\item{nu}{selection pressure. That means, number of offspring (lambda) is mue multiplied with nu. Default is \code{10}}
#' 			\item{dimension}{dimension number of the target function, default is \code{2}}
#' 			\item{mutation}{string of mutation type, default is \code{1}}
#' 			\item{sigmaInit}{initial sigma value (step size), default is \code{1.0}}
#' 			\item{nSigma}{number of different sigmas, default is \code{1}}
#' 			\item{tau0}{number, default is \code{0.0}. tau0 is the general multiplier.}
#' 			\item{tau}{number, learning parameter for self adaption, i.e. the local multiplier for step sizes (for each dimension).default is \code{1.0}}
#' 			\item{rho}{number of parents involved in the procreation of an offspring (mixing number), default is \code{"bi"}}
#' 			\item{sel}{number of selected individuals, default is \code{1}}
#' 			\item{stratReco}{Recombination operator for strategy variables. \code{1}: none. \code{2}: dominant/discrete (default). \code{3}: intermediate. \code{4}: variation of intermediate recombination. }
#' 			\item{objReco}{Recombination operator for object variables. \code{1}: none. \code{2}: dominant/discrete (default). \code{3}: intermediate. \code{4}: variation of intermediate recombination. }
#' 			\item{maxGen}{number of generations, stopping criterion, default is \code{Inf}}
#' 			\item{seed}{number, random seed, default is \code{1}}
#' 			\item{noise}{number, value of noise added to fitness values, default is \code{0.0}}
#' 			\item{lowerLimit}{number, lower limit for search space, default is \code{-1.0}}
#' 			\item{upperLimit}{number, upper limit for search space, default is \code{1.0}}
#' 			\item{verbosity}{defines output verbosity of the ES, default is \code{0}}
#' 			\item{plotResult}{boolean, specifies if results are plotted, default is \code{FALSE}}
#' 			\item{logPlotResult}{boolean, defines if plot results should be logarithmic, default is \code{FALSE}}
#' 			\item{term}{a string, defines which termination criterion should be used, default is \code{"iter"}}
#' 			\item{sigmaRestart}{number, value of sigma on restart, default is \code{0.1}}
#' 			\item{preScanMult}{initial population size is multiplied by this number for a pre-scan, default is \code{1}}
#' 			\item{globalOpt}{termination criterion on reaching a desired optimum value, default is \code{rep(0,dimension)}}
#'}
#'
#' @return This function returns a list with:\cr
#'	\code{par} parameters of the found solution\cr
#'	\code{value} target function value of the found solution\cr
#	\code{convergence} indicates successful completion when \code{0}\cr
#	\code{message}\cr
#	\code{hessian}\cr
#
#' @seealso \code{\link{optim}} \code{\link{spotOptim}}
#'
#' @export
###################################################################################################
spotOptimEs <- function(par #only for dimension TODO add par to start population
                      , fn
                      , gr=NULL
											, ...
                      , lower
                      , upper
                      , method=NULL
                      ,control=list()){					  
	if (length(par)==0) stop("dimension of par is null")
	con<-list(maxit=100, #CON: Internal List with defaults for control
		seed=1,
		mue = 10,
		nu = 10,
		mutation = 2,
		sigmaInit = 1.0,
		nSigma = 1,
		tau0 = 0.0,
		tau = 1.0,
		rho = "bi",
		sel = -1,
		stratReco = 1,
		objReco = 2,
		maxGen = Inf,
		noise = 0.0, 
		verbosity=0,
		plotResult=FALSE,
		logPlotResult=FALSE,
		term="iter",
		sigmaRestart = 0.1,
		preScanMult= 1);
	con[(namc <- names(control))] <- control
	control<-con;

	if(is.null(control$maxIter)){
		control$maxIter<-control$maxit
	}
	
	dimension=length(par)
	if(is.null(control$globalOpt)){
		control$globalOpt<-rep(0,dimension)
	}


	##### call SPOT	###############################			
	esResult<-spotAlgEs(mue=control$mue,
		nu=control$nu,
		dimension=dimension,
		mutation=control$mutation,
		sigmaInit=control$sigmaInit,
		nSigma=control$nSigma,
		tau0=control$tau0,
		tau=control$tau,
		rho=control$rho,
		sel=control$sel,
		stratReco=control$stratReco,
		objReco=control$objReco,
		maxGen=control$maxGen,
		maxIter=control$maxit,
		seed=control$seed,
		noise=control$noise,
		fName=fn,
		lowerLimit=lower,
		upperLimit=upper,
		verbosity=control$verbosity,
		plotResult=control$plotResult,
		logPlotResult=control$logPlotResult,
		term=control$term,
		sigmaRestart=control$sigmaRestart,
		preScanMult=control$preScanMult,
		globalOpt=control$globalOpt,		
		...)
	
	
	##### prepare Results	###############################		
	result <- list()
	result$par <- esResult$X
	result$value <- esResult$Y #TODO Y is real best, what about noisy value in case of noisy functions, or functions where noise is unknown?
	#result$counts <- #TODO log counts
	result$convergence <- 0 # 0 indicates successful completion	
	result$counts <- esResult$counts #number of function evaluations
	result
}

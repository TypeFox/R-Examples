#'	Latin Hypercube Sampling for Uncertainty and Sensitivity Analyses
#' 
#'	Generates the Latin Hypercube sampling for uncertainty and sensitivity analyses.
#' 
#'	A Latin Hypercube of size \code{N} is generated from the desired quantile distribution functions
#'
#'	The following methods are currently supported for generating the LHS: random LHS and 
#'	Huntington & Lyrintzis method for correcting the correlation matrix to be similar 
#'	to the prescribed by the option \code{COR} (see the arguments for description).
#'
#'	The specified \code{model} is run with the data from the LHS. If \code{repetitions}
#'	is set to more than one, the model will be run several times for each data point.
#'
#'	Partial rank correlation coefficients are estimated using code based on the \code{prcc}
#' function from the "sensitivity" package.
#'
#'	When the LHS function is called with no model (i.e., with argument
#'	\code{model=NULL}), it generates an incomplete object storing the Latin
#'	Hypercube samples, and allowing the user to run the simulation
#'	model independently. The method \code{\link{tell}} allows to pass the simulation
#'	results to the incomplete object.
#' 
#'  \code{tell} and \code{ask} are S3 generic methods for decoupling
#'  simulations and sensitivity measures estimations in the package
#'  `sensitivity'. In general, they are
#'  not used by the end-user for a simple \R model, but rather for an
#'  external computational code. The LHS object implements only the 
#'  \code{tell} method. For help on the other methods, see the help
#'  pages on the `sensitivity' package.
#'@param model The function to be run, representing the model or simulation.
#'		If \code{NULL}, no function is run and the object generated is incomplete, see also the \code{tell} method.
#'@param factors
#'		The names of the input variables (used for naming the 'data' data.frame and in plotting)
#'		Either a vector of strings or a single number representing the number of factors
#'@param N
#'		The size of the hypercube, i.e., how many samples are generated. Must be at least the number
#'		of factors plus 2.
#'@param q The quantile functions to be used. If only one is provided, it will be used for all parameters.
#'		Defaults to "qunif".
#'@param q.arg A list containing the arguments for the 'q' functions. Each parameter must be 
#' specified by a named list, containing all of the arguments for the quantile distribution. If unsupplied, 
#' default values	for the parameters are used.
#'@param res.names	Optional: what are the names of the model results? (Used mainly for plotting)
#'@param method Currently, two methods are supported. "random" generates a simple LH, with no modifications. 
#'		"HL" (the default) generates a random LH, and subsequently corrects the correlation matrix
#'		using the Huntington & Lyrintzis method.
#'@param opts
#'		Further options for the method used. The method HL supports the following options:
#'		`COR'	The desired correlation matrix between the model variables. If none is provided, the function will 
#'				generate a zero-correlation Latin Hypercube.
#'		`eps'   The tolerance between the prescribed correlation and the actual correlation present in the
#'				generated Latin Hypercube.
#'		`maxIt' The maximum number of iterations to be run for each factor. The default is set by a heuristic,
#'              but it might need some adjustments.
#'@param nboot Number of bootstrap replicates for calculating the PRCC.
#'@param repetitions The number of model repetitions to be run for a single data point. See the vignette on 
#'		stochastic models for details
#'@param cl	Cluster generated with the ``parallel'' library. May be of any type supported.
#'		If a cluster is provided, the model will be run in parallel or distributed across
#'		the cluster via clusterApply. No load balancing is provided, so the model results are
#'		reproducible.
#' 
#'		NOTE: You should manually export ALL objects required for the model to run, including the model
#'		function itself. See the help on \code{clusterExport} on package \code{parallel} for
#'		details.
#'@param x An LHS/PLUE object. For "tell", an incomplete LHS object (created with model=NULL)
#'@param y A data.frame containing the model responses
#'@param \dots Currently ignored
#'@section Warning:
#'	NOTE: the tell method from sensitivity objects (like 'fast99')
#'	modifies the object passed as argument as a side effect. This is 
#'	NOT the case with the LHS tell method.
#'@author Andre Chalom
#'@source Uses internal code originally published on package sensitivity, by Gilles Pujol, Bertrand Iooss, Alexandre Janon
#'@references
#'  McKay, M.D. and Beckman, R.J. 1979. A comparison of three methods for selecting 
#'  values of input variables in the analysis of output from a computer code, 
#'  \emph{Technometrics} 21: 239-244
#'
#'  Chalom, A. and Prado, P.I.K.L. 2012. Parameter space exploration of ecological models
#'  \emph{arXiv}:1210.6278 [q-bio.QM]
#'
#'@examples
#'completeLHS <- LHS(model=function(x) x[,1]+x[,2]*x[,3], factors=3, N=20)
#'incompleteLHS <- LHS(factors=5, N=30)
#'incompleteLHS <- tell(incompleteLHS, seq(1,30))
#'
#'\dontrun{
#'	new.cluster <- parallel::makePSOCKcluster(c("localhost", "localhost"))
#'	clusterLHS <- LHS(model=function(x) x[,1]/x[,2], factors=2, N=100, cl = new.cluster)
#'	stopCluster(new.cluster)
#'}
#' @export
LHS <-
	function (model=NULL, factors, N, q=NULL, q.arg=NULL, res.names=NULL, method=c("HL", "random"),
			  opts=list(), nboot=0, repetitions=1, cl = NULL) {
		# Input validation for common errors and "default" value handling:
		method = match.arg(method)
		my.opts = list(COR=0, eps=0.0005, maxIt=0)
		my.opts[names(opts)] <- opts
		if(is.numeric(factors) && length(factors) == 1) factors=paste("I", 1:factors, sep="")
		else if (!is.character(factors)) {
			stop("Error in function LHS: factors should be either a single number or a character vector")
		}
		if (!is.numeric(N) || length(N) != 1) {
			stop("Error in function LHS: N should be a single number");
		}
		if (N < length(factors) + 2) {
			stop("Error in function LHS: the number of points N must be at least the number of factors + 2");
		}
		# "Defaults" for q and q.arg
		if (is.null(q)) q=rep("qunif", length(factors)) 
		else if (length(q)==1)  q=rep(q, length(factors))
		if (is.null(q.arg)) q.arg =rep( list(list()), length(factors))
		else if (FALSE %in% sapply(q.arg, is.list)) q.arg <- rep(list(q.arg), length(factors))

		# Generates the hypercube data
		L <- as.data.frame(matrix(nrow=N, ncol=length(factors)));
		colnames(L) <- factors
		for (i in 1:length(factors)) 
			L[,i] <- sample(do.call(q[[i]], c(list(p = 1:N/N-1/N/2), q.arg[[i]])))
		# Corrects the correlation terms, for HL method for LHS
		if (method == "HL") {
			L <- LHScorcorr(L, COR = my.opts$COR, eps = my.opts$eps, maxIt = my.opts$maxIt); 
		}
		# Runs the actual model
		res <- internal.run(cl, model, L, repetitions)
		prcc <- internal.prcc(L, res, nboot)

		if (is.null(res.names) && ! is.na(res)) res.names <- paste("O", 1:dim(res)[2], sep="")
		X <- list(call=match.call(), N=N, data=L, factors=factors, q=q, q.arg=q.arg, 
				  opts = my.opts, model=model, res=res, prcc=prcc,
				  res.names=res.names);
		class(X) <- "LHS"
		return(X);
	}

internal.run <- function(cl, model, L, repetitions) {
	if (is.null(model)) {return(NA)}
	# First run, is independent of "repetitions"
	if (is.null(cl)) {
		tmp.res <- t(model(L));
		if(dim(tmp.res)[1] == 1) tmp.res = t(tmp.res)
	}
	else {
		tmp.res <- clusterRun(cl, model, L)
	}
	# and tells us the number of model outputs
	n.outs <- dim(tmp.res)[2]
	N <- dim(L)[1]
	res <- array(tmp.res, dim=c(N, n.outs, repetitions));
	if(repetitions> 1) for (i in 2:repetitions) {
		if (is.null(cl)) res[,,i] <- t(model(L))
		else res[,,i] <- clusterRun(cl, model, L)
	}
	return(res)
}

##Methods
#' @export
#' @rdname LHS
print.LHS <- function(x, ...) {
	  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
	  cat("Model:\n"); print (x$model);
	  cat("Factors:\n"); print (x$factors);
	  cat("Results:\n"); print (x$res.names);
	  cat("PRCC:\n"); print (x$prcc);
}

#' @export
#' @rdname LHS
tell <- function(x, y = NULL, ...)
	  UseMethod("tell")

#' @export
#' @rdname LHS
tell.LHS <- function (x, y, res.names=NULL, nboot=0, ...) {
	tmp.res <- t(y);
	if(dim(tmp.res)[1] == 1) tmp.res = t(tmp.res)
	# If this is the first "tell"
	if(all(is.na(get.results(x)))) {
		n.outs <- dim(tmp.res)[2]
		res <- array(tmp.res, dim=c(x$N, n.outs, 1));
	} else {
		this.repetition = get.repetitions(x)+1
		res <- array(get.results(x, get.mean=FALSE), dim=c(get.N(x), get.noutputs(x), this.repetition))
		res[,,this.repetition] = tmp.res
	}
	x$res <- res
	x$prcc <- internal.prcc(get.data(x), res, nboot)
	if (!is.null(res.names)) x$res.names <- res.names
	if (is.null(x$res.names)) x$res.names <- paste("O", 1:dim(res)[2], sep="")
	
	return(x)
}

internal.prcc <- function (L, res, nboot) {
	if (is.null(dim(res))) {return (NA)}
	# Reduces the res object to a 2-dimensional array
	res <- apply(res, c(1,2), mean)
	f <- function(r) pcc(L, r, nboot=nboot, rank=T)
	return(apply(res, 2, f))
}


#   Copyright (c) 2014-2015 by Martin Zaefferer, Cologne University of Applied Sciences

###################################################################################################
#' optimInterface
#'
#' This function is an interface fashioned like the \code{\link{optim}} function.
#' Unlike optim, it collects a set of bound-constrained optimization algorithms
#' with local as well as global approaches. It is used in the CEGO package
#' to solve the optimization problem that occurs during parameter estimation
#' in the Kriging model (based on Maximum Likelihood Estimation).
#' Note that this function is NOT applicable to combinatorial optimization problems. 
#'
#' The control list contains:\cr
#' \code{funEvals} stopping criterion, number of evaluations allowed for \code{fun}  (defaults to 100) \cr
#' \code{reltol} stopping criterion, relative tolerance  (default: 1e-6) \cr
#' \code{popsize} population size or number of particles  (default: \code{10*dimension}, where \code{dimension} is derived from the length of the vector \code{lower}). \cr
#' \code{restarts} the number of restarts to perform (Default: 0). Function evaluation budget will be split accordingly. Early convergence of optimization runs may lead to additional restarts. 
#' Violations of the provided budget may decrease the number of restarts.\cr
#' \code{method} will be used to choose the optimization method from the following list:\cr
#' "L-BFGS-B" - BFGS quasi-Newton: \code{stats} Package \code{optim} function\cr
#' "nlminb" - box-constrained optimization using PORT routines: \code{stats} Package \code{nlminb} function\cr
#' "DEoptim" - Differential Evolution implementation: \code{DEoptim} Package\cr
#' Additionally to the above methods, several methods from the package \code{nloptr} can be chosen. 
#' The complete list of suitable nlopt methods (non-gradient, bound constraints) is: \cr
#'	"NLOPT_GN_DIRECT","NLOPT_GN_DIRECT_L","NLOPT_GN_DIRECT_L_RAND",
#'	"NLOPT_GN_DIRECT_NOSCAL","NLOPT_GN_DIRECT_L_NOSCAL","NLOPT_GN_DIRECT_L_RAND_NOSCAL",
#'	"NLOPT_GN_ORIG_DIRECT","NLOPT_GN_ORIG_DIRECT_L","NLOPT_LN_PRAXIS",							
#'	"NLOPT_GN_CRS2_LM","NLOPT_LN_COBYLA",
#'	"NLOPT_LN_NELDERMEAD","NLOPT_LN_SBPLX","NLOPT_LN_BOBYQA","NLOPT_GN_ISRES"\cr\cr
#' All of the above methods use bound constraints.
#' For references and details on the specific methods, please check the documentation of the packages that provide them.\cr\cr
#' Furthermore, the user can choose to use a function instead of a string for the \code{method}.
#' The used function should have the same parameters and arguments as documented for this very function, i.e. \code{optimInterface}.
#' The \code{method} parameters for the call to the supplied function will be extracted from \code{control$method}.
#'
#' @param x is a point (vector) in the decision space of \code{fun}
#' @param fun is the target function of type \code{y = f(x, ...)}
#' @param lower is a vector that defines the lower boundary of search space
#' @param upper is a vector that defines the upper boundary of search space
#' @param control is a list of additional settings. See details.
#' @param ... additional parameters to be passed on to \code{fun}
#'
#' @return This function returns a list with:\cr
#'	\code{xbest} parameters of the found solution\cr
#'	\code{ybest} target function value of the found solution\cr
#'	\code{count} number of evaluations of \code{fun}
#'
#' @export
###################################################################################################
optimInterface<-function(x,fun,lower=-Inf,upper=Inf,control=list(),...){
	con<-list(funEvals=100 #CON: Internal List with defaults for control
			,method="L-BFGS-B"
      ,reltol=1e-6
			,popsize=NULL			
			,ineq_constr=NULL
			,verbosity=0
			,restarts=0)
	con[names(control)] <- control;
	control<-con;
	
	
	#INITIALIZE
	dim <- length(lower)
	if(is.null(control$popsize))
		control$popsize <- dim * 10
	budget <- control$funEvals / (control$restarts+1)
	sumevals <- 0	
	ymin <- Inf
	run <- TRUE
	method <- control$method
	
	if(!is.null(control$ineq_constr) &  !(method=="NLOPT_GN_ORIG_DIRECT" | method=="NLOPT_LN_COBYLA" ))
		warning("Constraint function passed to optimInterface. This is not supported with the chosen method.")
	
	if(length(x)==0) x <- runif(length(lower))*(upper-lower)+lower
	
	#LOOP OVER RESTARTS
	while(run){
		if (method=="L-BFGS-B"){
			res <- optim(par=x, fn=fun, method=method,lower=lower,upper=upper,control=list(maxit=budget,trace=control$verbosity),...)
			resval <- res$value
			respar <- res$par
			resevals <- res$counts[[1]] +res$counts[[1]] * 2 * dim
		}else if (method=="nlminb"){
			res <- nlminb(start=x, objective=fun, gradient=NULL, hessian=NULL, control=list(eval.max=budget,iter.max=budget,rel.tol=control$reltol,trace=control$verbosity),lower=lower,upper=upper,...)
			resval <- res$objective
			respar <- res$par
			resevals <- sum(res$evaluations)
		}else if (method=="DEoptim"){
			res <- DEoptim::DEoptim(fn=fun ,lower=lower,upper=upper,control=DEoptim.control(NP=control$popsize,itermax=floor((budget-control$popsize)/control$popsize),reltol=control$reltol,trace=FALSE),...)
			resval <- res$optim$bestval
			respar <- res$optim$bestmem
			resevals <- budget
		}else if (any(method==c("NLOPT_GN_DIRECT","NLOPT_GN_DIRECT_L","NLOPT_GN_DIRECT_L_RAND",
								"NLOPT_GN_DIRECT_NOSCAL","NLOPT_GN_DIRECT_L_NOSCAL","NLOPT_GN_DIRECT_L_RAND_NOSCAL",
								"NLOPT_GN_ORIG_DIRECT","NLOPT_GN_ORIG_DIRECT_L","NLOPT_LN_PRAXIS",							
								"NLOPT_GN_CRS2_LM","NLOPT_LN_COBYLA",
								"NLOPT_LN_NELDERMEAD","NLOPT_LN_SBPLX","NLOPT_LN_BOBYQA","NLOPT_GN_ISRES"))){ 
			opts=list(algorithm=method,maxeval=budget, ftol_rel=control$reltol, xtol_rel=-Inf, print_level=control$verbosity)	
			res <- nloptr::nloptr(x,fun,lb = lower,ub = upper, eval_g_ineq=control$ineq_constr,opts = opts,...)
			resval <- res$objective
			respar <- res$solution	
			resevals <- res$iterations
		}else{
			stop("The chosen optimization method used in optimInterface does not exist.")
		}
		
		#CONTROL RESTARTS
		sumevals <- sumevals+resevals
		if(resval < ymin){
			resultpar <- respar
			ymin <- resval
		}
		x <- lower+(upper-lower)*runif(dim)
		#Stop while loop when limit reached
		if(sumevals>=control$funEvals){
			run=FALSE
		}
		
	}
	result <- list()
	result$xbest <- resultpar
	result$ybest <- ymin
	result$count <- sumevals
	result

}

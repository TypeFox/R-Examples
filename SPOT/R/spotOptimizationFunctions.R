###################################################################################################
#' spotOptimizationInterface
#'
#' This function is an interface fashioned like the \code{\link{optim}} function. It is used in SPOT to access several different optimization methods.
#'
#' The control list contains:\cr
#' \code{fevals} stopping criterion, number of evaluations allowed for \code{fn}  (defaults to 100) \cr
#' \code{ineq_constr} defaults to NULL, else can be a function for an inequality constraint that is passed to nloptr\cr
#' \code{reltol} stopping criterion, relative tolerance  (defaults to 1e-6) \cr
#' \code{abstol} stopping criterion, absolute tolerance  (defaults to 1e-6) \cr
#' \code{popsize} population size or number of particles  (default depends on method) \cr
#' \code{restarts} whether or not to do restarts, default is FALSE   \cr
#' \code{vectorized} whether or not \code{fn} can evaluate multiple points at once, defaults to FALSE(only relevant for cmaes and pso methods) \cr
#' Please note that all settings will have to be passed to the actual optimization method. Not all of those
#' make use of the items listed above\cr\cr
#' Also note that the parameter \code{method} will be used to choose the optimization method from the following list:\cr
#' "lhs" - Latin Hypercube Sampling\cr
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
#' All of the above methods use bound constraints.
#' For references and details on the specific methods, please check the documentation of the packages that provide them.\cr\cr
#' Furthermore, the user can choose to use a function instead of a string for the \code{method}.
#' The used function should have the same parameters and arguments as documented for this very function, i.e. \code{spotOptimizationInterface}.
#'
#' @param par is a point (vector) in the decision space of \code{fn}
#' @param fn is the target function of type \code{y = f(x, ...)}
#' @param gr gradient function, not implemented yet
#' @param lower is a vector that defines the lower boundary of search space
#' @param upper is a vector that defines the upper boundary of search space
#' @param method is a string that describes which method is to be used, as implemented in this function. Else it can be a function, which is a custom
#' optimization function created by the user. See details.
#' @param control is a list of additional settings. See details.
#' @param ... additional parameters to be passed on to \code{fn}
#'
#' @return This function returns a list with:\cr
#'	\code{par} parameters of the found solution\cr
#'	\code{value} target function value of the found solution\cr
#'	\code{counts} number of evaluations of \code{fn}
#
#' @seealso \code{\link{spotOptim}}  \code{\link{spotOptimizationInterfaceMco}}
#'
#' @export
###################################################################################################
# TODO: par is a vector. should also be matrix in case of start population, or for restarts?
# TODO: this should also contain an LHS approach, to make the "largeDesign" stuff in the sequential step superfluous
spotOptimizationInterface<-function(par,fn,gr=NULL,lower,upper,method,control,...){
	#if (length(par)==0) stop("dimension of par is null")
	con<-list(fevals=100 #CON: Internal List with defaults for control
            ,abstol=1e-12
			,reltol=1e-12
			,popsize=20
			,ineq_constr=NULL
			,verbosity=0
			,restarts=FALSE
			,vectorized=FALSE)
	con[(namc <- names(control))] <- control;
	control<-con;
	
	#INITIALIZE
	dim = length(lower)
	sumevals <- 0	
	ymin <- Inf
	run=TRUE
	
	if(!is.null(control$ineq_constr) &  !(method=="NLOPT_GN_ORIG_DIRECT" | method=="NLOPT_LN_COBYLA" ))
		warning("Constraint function passed to spotOptimizationInterface. This is not supported with the chosen method.")
	
	# TODO infeasible par (start point)?
	
	if(length(par)==0) par <- runif(length(lower))*(upper-lower)+lower #TODO: in case of ineq constraint, only generate feasible startpoint
	
	#LOOP OVER RESTARTS
	while(run){
		if(is.function(method)){ #TODO test
			res <- method(par=par, fn=fn, gr=gr, lower=lower,upper=upper,control=control,...)
			resval <- res$value
			respar <- res$par
			resevals <- res$counts[[1]]
		}else if (method=="lhs"){
			res <- spotOptimLHS(par=par, fn=fn, lower=lower,upper=upper,control=list(fevals=control$fevals),...)
			resval <- res$value
			respar <- res$par
			#TODO: note: the following is only relevant for gradient free optimization
			resevals <- res$counts[[1]] 
		}else if (method=="optim-L-BFGS-B"){
			res <- optim(par=par, fn=fn, method="L-BFGS-B",lower=lower,upper=upper,control=list(maxit=control$fevals),...)
			resval <- res$value
			respar <- res$par
			#TODO: note: the following is only relevant for gradient free optimization
			resevals <- res$counts[[1]] +res$counts[[1]] * 2 * dim
		}else if (method=="BBoptim"){ 
			spotInstAndLoadPackages("BB")
			res <- BB::BBoptim(par=par, fn=fn,lower=lower,upper=upper,control=list(maxit=control$fevals,trace=FALSE),quiet=TRUE,...)
			resval <- res$value
			respar<- res$par
			resevals <- res$feval
		}else if (method=="pso"){ #TODO note: pso would be able to do restarts internally, integrate?
			spotInstAndLoadPackages("pso")
			res <- pso::psoptim(par=par, fn=fn,lower=lower,upper=upper,control=list(vectorize=control$vectorized, maxf=control$fevals,reltol=control$reltol),...)
			resval <- res$value
			respar <- res$par
			resevals <- res$counts[[1]]	
		}else if (method=="cmaes"){ 
			spotInstAndLoadPackages("cmaes");
			if(control$vectorized)fun=function(x)fn(t(x))   #cmaes has columns for observations, rows for parameters
			else fun <- fn
			res <- cmaes::cma_es(par, fun,lower=lower,upper=upper,control=list(vectorized=control$vectorized,maxit=ceiling(control$fevals / (4 + floor(3 * log(dim))))),...) 
			if(is.null(res$par)){#error handling, see bug example in start.runit.test at end of file, reason unclear
				resval <- fn(par,...)
				respar <- par
				resevals <- 1	
			}
			else{
				resval <- res$value
				respar <- res$par
				resevals <- res$counts[[1]]	
			}			
		}else if (method=="genoud"){ 
			spotInstAndLoadPackages("rgenoud")
			res <- rgenoud::genoud(nvars=dim,starting.values=par,fn<-fn,boundary.enforcement=1,Domains=cbind(lower,upper),print.level=0,pop.size=control$popsize,max.generations=floor(control$fevals/control$popsize),wait.generations=floor(control$fevals/control$popsize),...)
			resval <- res$value
			respar <- res$par
			resevals <- res$generations * control$popsize		
		}else if (method=="DEoptim"){
			spotInstAndLoadPackages("DEoptim")
			res <- DEoptim::DEoptim(fn=fn ,lower=lower,upper=upper,control=DEoptim::DEoptim.control(NP=control$popsize,itermax=floor((control$fevals-control$popsize)/control$popsize),reltol=control$reltol,trace=FALSE),...)
			resval <- res$optim$bestval
			respar <- res$optim$bestmem
			resevals <- control$fevals
		}else if (method=="bobyqa"){ 
			spotInstAndLoadPackages("minqa")
			res <- minqa::bobyqa(par,lower=lower,upper=upper,fn ,control=list(maxfun=control$fevals,iprint=FALSE),...)
			resval <- res$fval
			respar <- res$par
			resevals <- res$feval
		}else if (method=="hjkb"){ 
			spotInstAndLoadPackages("dfoptim")
			res <- dfoptim::hjkb(par=par,fn=fn, lower=lower, upper=upper ,control=list(maxfeval=control$fevals),...)
			resval <- res$value
			respar <- res$par
			resevals <- res$feval	
		}else if (method=="GenSA"){ 
			spotInstAndLoadPackages("GenSA");
			res <- GenSA::GenSA(par=,fn=fn, lower=lower, upper=upper ,control=list(max.call=control$fevals),...)
			resval <- res$value
			respar <- res$par
			resevals <- res$counts[[1]]
		}else if (any(method==c("NLOPT_GN_DIRECT","NLOPT_GN_DIRECT_L","NLOPT_GN_DIRECT_L_RAND",
								"NLOPT_GN_DIRECT_NOSCAL","NLOPT_GN_DIRECT_L_NOSCAL","NLOPT_GN_DIRECT_L_RAND_NOSCAL",
								"NLOPT_GN_ORIG_DIRECT","NLOPT_GN_ORIG_DIRECT_L","NLOPT_LN_PRAXIS",							
								"NLOPT_GN_CRS2_LM","NLOPT_LN_COBYLA","NLOPT_LN_NEWUOA_BOUND",
								"NLOPT_LN_NELDERMEAD","NLOPT_LN_SBPLX","NLOPT_LN_BOBYQA","NLOPT_GN_ISRES"))){ 
			opts=list(algorithm=method,maxeval=control$fevals, ftol_rel=control$reltol, xtol_rel=-Inf)	
			res <- nloptr::nloptr(par,fn,lb = lower,ub = upper, eval_g_ineq=control$ineq_constr,opts = opts,...)
			resval <- res$objective
			respar <- res$solution	
			resevals <- res$iterations
			# TODO: converges not to feasible solutions ?
		}else{
			stop("The chosen optimization method used in spotOptimizationInterface does not exist")
		}
		#CONTROL RESTARTS
		if(control$restarts){
			sumevals <- sumevals+resevals
			if(resval < ymin){
				resultpar <- respar
				ymin <- resval
			}
			par <- lower+(upper-lower)*runif(dim)
			#Stop while loop when limit reached
			if(sumevals>=control$fevals){
				run=FALSE
			}
		}else{ #IN CASE OF NO RESTARTS
			ymin <- resval
			run <- FALSE
			resultpar <- respar
			sumevals <- resevals			
		}
	}
	result <- list()
	result$par <- resultpar
	result$value <- ymin
	result$counts <- sumevals
	result

}

###################################################################################################
#' spotOptimizationInterfaceMco
#'
#' This function is an interface fashioned like the \code{\link{optim}} function. It is used in SPOT to access several different multi-criteria optimization methods.
#'
#' The control list contains:\cr
#' \code{fevals} stopping criterion, number of evaluations allowed for \code{fn}  (defaults to 100) \cr
# \code{reltol} stopping criterion, relative tolerance  (defaults to 1e-6) \cr
# \code{abstol} stopping criterion, absolute tolerance  (defaults to 1e-6) \cr
#' \code{popsize} population size or number of particles  (default depends on method) \cr
#' \code{restarts} whether or not to do restarts, default is FALSE   \cr
# \code{vectorized} whether or not \code{fn} can evaluate multiple points at once, defaults to FALSE (only relevant for cmaes and pso methods) \cr
#' Also note that the parameter \code{method} will be used to choose the optimization method from the following list:\cr
#' "nsga2" - the nsga2 function from the \code{mco} Package \cr
#' "sms-emoa" - The basic sms-emoa in the \code{SPOT} package\cr
#' All of the above methods use bound constraints.
#' For references and details on the specific methods, please check the documentation of the packages that provide them.\cr\cr
#' Furthermore, the user can choose to use a function instead of a string for the \code{method}.
#' The used function should have the same parameters and arguments as documented for this very function, i.e. \code{spotOptimizationInterfaceMco}.
#'
#' @param par is a point (vector) in the decision space of \code{fn}, par is not used (yet)
#' @param fn is the target function of type \code{y = f(x, ...)}
#' @param gr gradient function, gr is not used (yet)
#' @param lower is a vector that defines the lower boundary of search space
#' @param upper is a vector that defines the upper boundary of search space
#' @param method is a string that describes which method is to be used, as implemented in this function. Else it can be a function, which is a custom
#' optimization function created by the user. See details.
#' @param control is a list of additional settings. See details.
#' @param ref reference point. Please provide this even with methods that do not use it (e.g. "nsga2"), to specify the dimension of the objective space.
#' @param ... additional parameters to be passed on to \code{fn}
#'
#' @return This function returns a list with:\cr
#'	\code{par} parameters of the found solutions, e.g. the Pareto set\cr
#'	\code{value} target function values of the found solutions, e.g. the Pareto front\cr
#'	\code{counts} number of evaluations of \code{fn}
#
#' @seealso \code{\link{spotOptimizationInterface}} \code{\link{spotOptim}}
#'
#' @export
###################################################################################################
# TODO: par is a vector. should also be matrix in case of start population, or for restarts?
# TODO: this should also contain an LHS approach, to make the "largeDesign" stuff in the sequential step superfluous
spotOptimizationInterfaceMco<-function(par,fn,gr=NULL,lower,upper,method,control,ref,...){
	#if (length(par)==0) stop("dimension of par is null")
	con<-list(fevals=100 #CON: Internal List with defaults for control
            ,abstol=1e-12
			,reltol=1e-12
			,popsize=20
			,verbosity=0
			,sbx.n=15, sbx.p=0.7,pm.n=25, pm.p=0.3
			,restarts=FALSE
			,vectorized=FALSE)
	con[(namc <- names(control))] <- control;
	control<-con;
	
	#INITIALIZE
	dimi = length(lower)
	dimo = length(ref)
	sumevals <- 0	
	set<-NULL
	front<-NULL
	run=TRUE
	#LOOP OVER RESTARTS
	while(run){
		if(is.function(method)){ #TODO untested
			res <- method(par=par, fn=fn, gr=gr, lower=lower,upper=upper,control=control,...)
			resval <- res$value
			respar <- res$par
			resevals <- res$counts[[1]]
		}else if (method=="nsga2"){
			spotInstAndLoadPackages("mco");		
			psize = control$popsize
			if((psize%%4)!=0)psize=psize+4-(psize%%4) #make sure that psize is multiple of 4 for nsga2
			r1 <- nsga2(fn, dimi, dimo,
			   generations=floor(control$fevals/psize), popsize=psize,
			   lower.bounds=lower, 
			   upper.bounds=upper, 
			   cprob=control$sbx.p,
			   cdist=control$sbx.n,
			   mprob=control$pm.p,
			   mdist=control$pm.n,  ...)
			resval <- r1$value;		
			respar <- r1$par;
			resevals <- control$fevals	
		}
		else if (method=="sms-emoa"){
			psize = control$popsize
			fun <- function(x,...){if(any(is.na(x))){rep(NA,dimo)}else{fn(x)}} #wrapper. sms-emoa uses na values to check o-dim
			r1 <- spotSmsEmoa(fun,
			   lower=lower,
			   upper=upper,
			   control=list(mu=psize,sbx.n=control$sbx.n,sbx.p=control$sbx.p,pm.n=control$pm.n,pm.p=control$pm.p,maxeval=control$fevals),...)
			#ndom <- !is_dominated(r1$Y)
			#resval <- t(r1$Y[,ndom])
			#respar <- t(r1$X[,ndom])
			resval <- t(r1$value);		
			respar <- t(r1$par);
			resevals <- control$fevals
		}
		else{
			stop("The chosen optimization method used in spotOptimizationInterfaceMco does not exist")
		}
		#CONTROL RESTARTS
		if(control$restarts){
			sumevals <- sumevals+resevals
			set <- rbind(set,respar)
			front <- rbind(front,resval)
			#Stop while loop when limit reached
			if(sumevals>=control$fevals){
				run=FALSE
			}
		}else{ #IN CASE OF NO RESTARTS
			run <- FALSE
			set <- respar
			front <- resval
			sumevals <- resevals			
		}
	}
	result <- list()
	result$par <- set
	result$value <- front
	result$counts <- sumevals
	result
}

###################################################################################################
#' spotOptimLHS
#'
#' This function is an interface to Latin Hypercube Sampling (LHS) fashioned like the \code{\link{optim}} function. 
#' That means, LHS is performed to optimize a target function, i.e. returning the sample with the 
#' lowest function value. 
#'
#' The control list contains:\cr
#' \code{fevals} number of design points created \cr
#' \code{retries} number of designs created during creation of a well spread design\cr
#' \code{vectorized} whether or not \code{fn} can evaluate multiple points at once , defaults to FALSE
#'
#' @param par is a point (vector) in the decision space of \code{fn}. Points in par will be added to the design created by LHS.
#' @param fn is the target function of type \code{y = f(x, ...)}
#' @param gr gradient function, gr is not used (yet)
#' @param lower is a vector that defines the lower boundary of search space
#' @param upper is a vector that defines the upper boundary of search space
#' @param control is a list of additional settings. See details.
#' @param ... additional parameters to be passed on to \code{fn}
#'
#' @return This function returns a list with:\cr
#'	\code{par} parameters of the found solutions, e.g. the Pareto set\cr
#'	\code{value} target function values of the found solutions, e.g. the Pareto front\cr
#'	\code{counts} number of evaluations of \code{fn}
#
#' @seealso \code{\link{spotOptimizationInterface}} \code{\link{spotOptim}}
#'
#' @export
###################################################################################################
spotOptimLHS<-function(par,fn,gr=NULL,lower,upper,control,...){
	#if (length(par)==0) stop("dimension of par is null")
	con<-list(fevals=1000 #CON: Internal List with defaults for control
			,retries=50
			,vectorized=FALSE)
	con[(namc <- names(control))] <- control;
	control<-con;
	
	npar <- nrow(par)
	npar <- ifelse(is.null(npar),1,npar)
		
	size <- control$fevals - npar
	dimension <- length(lower)
	best <- spotNormDesign(dimension,size,calcMinDistance=control$retries>1,nested=par);
	
	if (control$retries>1) {
		for (i in 1:(control$retries-1)) {
			tmpDes <- spotNormDesign(dimension,size,calcMinDistance=TRUE,nested=par);
			## maximize minimal distance
			if (tmpDes$minDistance > best$minDistance)
				best <- tmpDes;
		}
	}
	#scale design into bounds
	best<- t(lower+ t(best$design) * (upper-lower))
	#par is added to design
	par<-rbind(best,par)

	#evaluate design
	if(control$vectorized){
		value <- fn(par,...)
	}else{
		value <- apply(par,1,fn,...)
	}
	#determine best point
	ind <- which.min(value)
	#return!
	result <- list()
	result$par <- par[ind,]
	result$value <- value[ind]
	result$counts <- length(value)
	result
}


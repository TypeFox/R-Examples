turboem <- function(par, fixptfn, objfn = NULL, method = c("em", "squarem", "pem", "decme", "qn"), boundary = NULL, pconstr = NULL, project = NULL, parallel = FALSE, ..., control.method = replicate(length(method),list()), control.run = list()) {
	# method = vector of one or more methods
	# control.method = list containing control parameters if length(method)==1
	#				 = list(list1, list2, ...) when method = c(method1, method2, ...), where list1 is the list of control parameters for method1, list2 is the list of control parameters for method2, ...
	# control.run = list of control parameters for convergence and stopping rule across each method
	
	## allow partial matching of method and all lower case
	method <- tolower(method)
	method <- match.arg(method, several.ok=TRUE)
	
	if(class(control.method) != "list")
		stop("'control.method' must be of class 'list'")
	if(length(method)==1 & class(control.method[[1]]) != "list")
		control.method <- list(control.method)
	if(length(control.method) != length(method))
		stop("The number of components of 'control.method' must be equal to the number of methods being compared")
	for(j in seq_along(method))
		if(class(control.method[[j]]) != "list")
			stop("each component of 'control.method' must be of class 'list'")
			
	if(missing(objfn)) objfn <- NULL
			
	control.run.default <- list(convtype = "parameter", tol = 1e-07, stoptype = "maxiter", maxiter = 1500, maxtime = 60, convfn.user = NULL, stopfn.user = NULL, trace = FALSE, keep.objfval = FALSE)
	namc <- names(control.run)
	if(!all(namc %in% names(control.run.default))) {stop("unknown names in control.run: ", namc[!(namc %in% names(control.run.default))])}
	control.run <- modifyList(control.run.default, control.run)
	
	pars <- matrix(NA, length(method), length(par))
	rownames(pars) <- method
	objfeval <- fpeval <- itr <- value.objfn <- vector("numeric", length(method))
	convergence <- fail <- vector("logical", length(method))
	errors <- vector("character", length(method))
	runtime <- matrix(NA, length(method), 3)
	colnames(runtime) <- c("user", "system", "elapsed")
	if(!control.run$keep.objfval) {
		trace.objfval <- NULL
	} else {
		trace.objfval <- vector("list", length(method))
	}
	# Modified (JFB 30Jan2012)
	if(parallel) {
		resall <- foreach(j = seq_along(method), .packages="turboEM", .errorhandling = "pass") %dopar% accelerate(par=par, fixptfn=fixptfn, objfn=objfn, boundary=boundary, method=method[j], pconstr=pconstr, project=project, ..., control=modifyList(control.run, control.method[[j]]))
	} else {
		resall <- foreach(j = seq_along(method), .packages="turboEM", .errorhandling = "pass") %do% accelerate(par=par, fixptfn=fixptfn, objfn=objfn, boundary=boundary, method=method[j], pconstr=pconstr, project=project, ..., control=modifyList(control.run, control.method[[j]]))		
	}
	for(j in seq_along(method)) {
		if(control.run$trace) cat(paste("Method ", j, ":  \n", sep=""))
		res <- resall[[j]]
		if(!"turboem" %in% class(res)) {
			fail[j] <- TRUE
			##errors[j] <- res[[1]]
			errors[j] <- as.character(res)
			pars[j,] <- NA
			itr[j] <- NA
			fpeval[j] <- NA
			objfeval[j] <- NA
			value.objfn[j] <- NA
			convergence[j] <- NA
			if(control.run$keep.objfval) {
				trace.objfval[[j]] <- NA
			}
		} else {
			errors[j] <- NA
			pars[j,] <- res$par
			itr[j] <- res$itr
			fpeval[j] <- res$fpeval
			objfeval[j] <- res$objfeval
			value.objfn[j] <- res$value.objfn
			convergence[j] <- res$convergence
			runtime[j,] <- round(res$runtime[1:3], 3)
			if(control.run$keep.objfval) {
				trace.objfval[[j]] <- res$trace.objfval
			}
			control.method[[j]] <- res$control
		}
	}
	lst <- list(...)
	arg.names <- names(lst)
	if(is.null(objfn)) {
		obfjn.return <- NULL
	} else {
		obfjn.return <- objfn
		formals(obfjn.return)[arg.names] <- lst
	}
	formals(fixptfn)[arg.names] <- lst
	if(!is.null(pconstr)) {
		sel <- names(formals(pconstr))[names(formals(pconstr)) %in% arg.names]
		formals(pconstr)[sel] <- lst[sel]
	}
	ret <- list(fail=fail, value.objfn=value.objfn, itr=itr, fpeval=fpeval, objfeval=objfeval, convergence=convergence, runtime=runtime, errors=errors, pars=pars, method=method, trace.objfval=trace.objfval, control.method=control.method, control.run=control.run, fixptfn=fixptfn, objfn=obfjn.return, pconstr=pconstr, project=project)
	class(ret) <- "turbo"
	return(ret)
}
########################################
########################################
## Function to implement a particular acceleration scheme
########################################
########################################

accelerate <- function(par, fixptfn, objfn, method = c("em", "squarem", "pem", "decme", "qn"), boundary = NULL, pconstr = NULL, project = NULL, ..., control = list()) {
	# par = starting value of parameter vector
	# fixptfn = fixed-point iteration F(x)
	# for which the solution: F(x*) = x* is sought
	# objfn = underlying objective function which is minimized at x*
	# boundary = NULL unless method = "decme"
	# pconstr = NULL unless there are constraints on par
	# control = list of general control parameters for convergence and stopping criteria, as well as algorithm-specific control parameters
	
	env <- new.env()
	control.default <- list(convtype = "parameter", tol = 1e-07, stoptype = "maxiter", maxiter = 1500, maxtime = 60, convfn.user = NULL, stopfn.user = NULL, trace = FALSE, keep.objfval = FALSE, use.pconstr = TRUE)
	## added use.pconstr (8Feb2012 JFB): this allows user to specify that pconstr not be used even if it is included (useful for benchmark studies where you want to compare including it to not including it)
	
	## allow partial matching of method and all lower case
	method <- tolower(method)
	method <- match.arg(method)
	
	## if par is a matrix, convert to a vector
	if(is.matrix(par)) par <- as.vector(par)	
	
	## control parameters for each algorithm
	if(method=="em") {
		ctrl <- control.default
	} else if(method=="pem") {
		control.parabolicEM <- list(l=10, h=0.1, a=1.5, version="geometric")
		ctrl <- c(control.default, control.parabolicEM)
		if(is.null(objfn)) {stop("objfn required for method = 'pem' \n")}
	} else if(method=="decme") {
		if(is.null(boundary)) {stop("decme requires boundary function \n")}
		if(is.null(objfn)) {stop("objfn required for method = 'decme' \n")}
		control.dynamicECME <- list(version = "2s", tol_op = 0.01)
		ctrl <- c(control.default, control.dynamicECME)
	} else if(method=='qn'){
		if(is.null(objfn)) {stop("objfn required for method = 'qn' \n")}
		if(is.null(pconstr)) {
			warning("If the parameter space is constrained, then pconstr should be provided for method='qn'")
		}
		control.quasiNewton <- list(qn = 2)
		ctrl <- c(control.default, control.quasiNewton)
	} else if(method=="squarem") {
		control.squarem <- list(K = 1, square = TRUE, version = 3, step.min0 = 1, step.max0 = 1, mstep = 4, kr = 1, objfn.inc = 1, trace = FALSE)	
		ctrl <- c(control.default, control.squarem)
	}
	namc <- names(control)
	if(!all(namc %in% names(ctrl))) {stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])}
	ctrl <- modifyList(ctrl, control)
	
	convtype <- ctrl$convtype
	tol <- ctrl$tol
	stoptype <- ctrl$stoptype
	maxiter <- ctrl$maxiter
	maxtime <- ctrl$maxtime
	convfn.user <- ctrl$convfn.user
	stopfn.user <- ctrl$stopfn.user
	trace <- ctrl$trace
	keep.objfval <- ctrl$keep.objfval
	
	## to store objective function value at each iteration
	trace.objfval <- NULL
	
	## user-defined constraints on parameters
	# modified 8Feb2012 (JFB)
	if(is.null(pconstr) | !ctrl$use.pconstr) {
	# if(is.null(pconstr)) {
		pconstr <- function(par) TRUE
	}
	
	## default projection function: no projection is executed 
	if(is.null(project))
		project <- function(par) par
	environment(project) <- env
	
	## convergence criterion
	if(is.null(convfn.user)) {
		if(convtype=="parameter") {
			# Modified (JFB 30Jan2012)
			convfn <- function(old, new) {sqrt(crossprod(new-old)) < tol}
			# convfn <- function(old, new) {
				# diff <- as.numeric(sqrt(crossprod(new-old)))
				# isTRUE(all.equal(diff, tol, tolerance = tol))
			# }
		} else if(convtype=="objfn") {
			# Modified (JFB 30Jan2012)
			convfn <- function(old, new) {abs(new-old) < tol}
			# convfn <- function(old, new) {
				# diff <- max(abs(new-old))
				# isTRUE(all.equal(diff, tol, tolerance = tol))
			# }
		}
	}
	if(!is.null(convfn.user)) {
		convfn <- convfn.user
	}
	if(!convtype %in% c('parameter', 'objfn')) {
		stop("\n convtype must be one of c('parameter', 'objfn') \n")
	}
	if(convtype=="objfn" & is.null(objfn)) {
		stop("\n objfn required if convtype='objfn' \n")
	}
	environment(convfn) <- env
	
	## criterion for upper bound stopping the algorithm (e.g. maximum number of iterations or maximum runtime)
	if(stoptype=="maxiter" & is.null(stopfn.user)) {
		stopfn <- function() {iter == maxiter}
	} else if(stoptype=="maxtime" & is.null(stopfn.user)) {
		stopfn <- function() {elapsed.time >= maxtime}
	} else if(!is.null(stopfn.user)) {
		ctrl$stoptype <- "user"
		stopfn <- stopfn.user
		environment(stopfn) <- env
	} 
	
	convergence <- FALSE
	iter <- 0
	fpeval <- 0
	objfeval <- 0
	STOP <- FALSE
	start.time <- proc.time()
	
	########################################
	########################################
	## initialize algorithm ################
	########################################
	########################################
	if(method=="em") {
		input <- list(par=par, objfval=NULL)
		mainfun <- bodyEM
		rm(par)
	} else if(method=="pem") {
		## run EM algorithm for l iterations to get starting values
		for(i in 1:ctrl$l) {
			par <- fixptfn(par, ...)
			fpeval <- fpeval + 1
			if(i==ctrl$l-2) {
				p0 <- par
			} else if(i==ctrl$l-1) {
				p1 <- par
			}
		}
		p2 <- par
		input <- list(p0=p0, p1=p1, par=p2, objfval=objfn(p2, ...))
		objfeval <- objfeval + 1
		ctrl.alg <- ctrl[c("h", "a", "version")]
		mainfun <- bodyParaEM
		rm(p0, p1, p2, par)
	} else if(method=="decme") {
		if(ctrl$version=="2s") {
			objfval.old <- objfn(par, ...)
			objfeval <- objfeval + 1
		} else objfval.old <- 1e100
		input <- list(par=par, objfval.old=objfval.old, objfval=1e100)
		ctrl.alg <- c(ctrl[c("version", "tol_op")], dm=length(par))
		mainfun <- bodyDECME
		rm(par)
	} else if(method=='qn') {
		if(ctrl$qn > length(par)) {
			stop("For quasi-Newton, control parameter 'qn' must not be larger than the dimension of 'par'.")
		}
		U <- NULL
		for(i in 1:ctrl$qn)
		{
			pnew <- fixptfn(par, ...)
			fpeval <- fpeval + 1
			U    <- cbind(U,pnew-par)
			par  <- pnew
		}
		pnew  <- fixptfn(par, ...)
		fpeval <- fpeval + 1
		U     <- cbind(U[,-1],pnew-par)
		pnew2 <- fixptfn(pnew, ...)
		fpeval <- fpeval + 1
		V <- cbind(U[,-1],pnew2-pnew)
		
		input <- list(par=par, pnew=pnew, pnew2=pnew2, V=V, U=U)
		ctrl.alg <- ctrl["qn"]
		mainfun <- bodyQuasiNewton
		rm(par, pnew, pnew2, V, U)
	} else if(method=="squarem") {
		if (ctrl$K > 1 & !(ctrl$version %in% c("rre", "mpe"))) 
			ctrl$version <- "rre"
		if (ctrl$K == 1 & !(ctrl$version %in% c(1, 2, 3))) 
			ctrl$version <- 3
		if (!is.null(objfn)) {
			if (ctrl$K == 1) { 
				mainfun <- bodySquarem1
				if (trace) 
					cat("Squarem-1 \n")
				if (is.null(objfn)) 
					stop("\n squarem2 should be used if objective function is not available \n")
				input <- list(par=par, objfval=objfn(par, ...), step.max=ctrl$step.max0, step.min=ctrl$step.min0)
				objfeval <- objfeval + 1
				ctrl.alg <- ctrl[c("version", "step.max0", "mstep", "objfn.inc")]
				rm(par)
			} else if (ctrl$K > 1 | ctrl$version %in% c("rre", "mpe")) {
				mainfun <- bodyCyclem1
				if (trace) 
					cat("Cyclem-1 \n")
				if (is.null(objfn)) 
					stop("\n cyclem2 should be used if objective function is not available \n")
				if (ctrl$K > length(par)) {
					stop("K is too large.  Decrease it.  \n")
					return()
				}
				input <- list(par=fixptfn(par, ...), objfval=objfn(par, ...))
				fpeval <- fpeval + 1
				objfeval <- objfeval + 1
				ctrl.alg <- ctrl[c("version", "K", "square", "objfn.inc")]
				rm(par)
			}
		} else {
			if (ctrl$K == 1) {
				mainfun <- bodySquarem2
				if (trace) 
					cat("Squarem-2 \n")
				input <- list(par=par, objfval=NULL, step.max=ctrl$step.max0, step.min=ctrl$step.min0)
				ctrl.alg <- ctrl[c("version", "kr", "step.max0", "mstep")]
				rm(par)
			} else if (ctrl$K > 1 | ctrl$version %in% c("rre", "mpe")) {
				mainfun <- bodyCyclem2
				if (trace) 
					cat("Cyclem-2 \n")
				if (ctrl$K > length(par)) {
					cat("K is too large.  Decrease it.  \n")
					return()
				}
				input <- list(par=fixptfn(par, ...), objfval=NULL)
				fpeval <- fpeval + 1
				ctrl.alg <- ctrl[c("version", "K", "square", "kr")]
				rm(par)
			}
		}		
	}
	environment(mainfun) <- env
	
	########################################
	########################################
	## Set up and iterate ##################
	########################################
	########################################
	par <- input$par
	if(is.null(input$objfval)) {
		if(convtype=="objfn" | keep.objfval) {
			if(is.null(objfn)) stop("objfn must be supplied if keep.objfval=TRUE")
			objfval.old <- objfn(input$par, ...)
			objfeval <- objfeval + 1
		}
	} else {
		objfval.old <- input$objfval
	}
	## trace of objective function value
	if(keep.objfval) {
		##trace.objfval <- rbind(trace.objfval, c((proc.time()-start.time)[3], objfval.old))
		trace.objfval <- c(trace.objfval, objfval.old)
	}
	## iterate through algorithm
	iter.starttime <- proc.time()-start.time
	while(!STOP) {
		iter <- iter+1
		
		## main body of the algorithm
		itrout <- mainfun(input, ctrl=ctrl.alg, fixptfn=fixptfn, objfn=objfn, pconstr=pconstr, project=project, iter=iter, fpeval=fpeval, objfeval=objfeval, boundary=boundary, ...)
		output <- itrout$output
		p.new <- output$par
		objfval.new <- output$objfval
		fpeval <- itrout$fpeval
		objfeval <- itrout$objfeval
		
		if(trace) {
			cat("Objective fn: ", objfval.new, "  Iter: ", 
					round(iter), "  diff", sqrt(crossprod(p.new-par)), "\n")
		}
		
		## inputs necessary to check if convergence been achieved
		if(convtype=="parameter") {
			old <- par
			new <- p.new
		} else if(convtype=="objfn") {
			# if(is.null(objfval.old)) {
			# objfval.old <- objfn(par, ...)
			# objfeval <- objfeval + 1
			# }
			old <- objfval.old
			if(is.null(objfval.new)) {
				objfval.new <- objfn(p.new, ...)
				objfeval <- objfeval + 1
			} 
			new <- objfval.new
		}
		
		## history of objective function value, if desired
		if(keep.objfval) {
			if(is.null(objfval.new)) {
				objfval.new <- objfn(input$par, ...)
				objfeval <- objfeval + 1
			}
			##trace.objfval <- rbind(trace.objfval, c((proc.time()-start.time)[3], objfval.new))
			trace.objfval <- c(trace.objfval, objfval.new)
		}
		
		## has convergence been achieved?			
		if(convfn(old=old, new=new)) {
			convergence <- TRUE
			break
		}
		par <- p.new
		objfval.old <- objfval.new
		input <- output
		
		## should the algorithm be stopped due to other constraints (e.g. on run time or number of iterations)?
		elapsed.time <- (proc.time()-start.time)[3] ## elapsed time
		if(stopfn()) {
			STOP <- TRUE
		}
	}
	if(is.null(objfval.new) & !is.null(objfn)) {
		objfval.new <- objfn(p.new, ...)
		objfeval <- objfeval + 1
	} else if(is.null(objfval.new) & is.null(objfn)) {
		objfval.new <- NA
	}
	
	if(keep.objfval) {
		time.before.iter <- iter.starttime
		time.per.iter <- (proc.time() - start.time - iter.starttime)/iter
		trace.objfval <- list(time.before.iter=time.before.iter, time.per.iter=time.per.iter, trace=trace.objfval)
	}
	ret <- list(par=p.new, value.objfn=objfval.new, itr=iter, fpeval=fpeval, objfeval=objfeval, convergence=convergence, method=method, runtime=proc.time()-start.time, keep.objfval=keep.objfval, trace.objfval=trace.objfval, control=ctrl[!names(ctrl) %in% names(control.default)])
	class(ret) <- "turboem"
	return(ret)
}

########################################
########################################
## Additional functions needed for different algorithms
########################################
########################################

## function used within DECME-2
line.search <- function(pnew, p1, p0, llnew, tol_op, boundary, objfn, ...) {
	# maximize
	# pnew: current estimation
	# p1-p0 defines the searching direction
	dr <- p1-p0
	temp <- boundary(pnew, dr)
	optimiter <- 0
	fct <- function(rho, par, dr, ...) {
		# to count number of objective function evaluations
		optimiter <<- optimiter + 1
		return(objfn(par+rho*dr, ...))
	}
	res <- try(optimize(fct, temp, par=pnew, dr=dr, tol=tol_op, ...), silent=TRUE)
	if( class(res)!="try-error"){
		if(is.finite(res$objective)){
			if(res$objective < -llnew){
				pnew.tmp <- pnew + res$minimum*dr
				pnew <- pnew.tmp
				llnew <- -res$objective
			}
		}
	}
	return(list(pnew=pnew, llnew=llnew, optimiter=optimiter))
}
## functions used within DECME-2s
inv2 <- function(a) {
	return(matrix(c(a[2,2], -a[1,2], -a[2,1], a[1,1])/(a[1,1]*a[2,2]-a[1,2]*a[2,1]),2,2))
}
search <- function(pnew, dr, dm, boundary, objfn, ...) {
	bd <- boundary(pnew, dr)
	
	if(bd[1]==0) {
		alpha <- 0
		lltheta <- -Inf
		theta <- rep(0, dm)
		objfiter <- 0
	} else {
		alpha <- min(1, bd[2]*0.9)
		theta <- pnew + alpha*dr
		lltheta <- try(-objfn(theta, ...), silent=T)
		objfiter <- 1
	}
	
	return(list(res=c(alpha, lltheta, theta), dr=dr, objfiter=objfiter))
}
optimize.2d.app=function(pnew, pold, par, llnew, llold, llp, dm, boundary, objfn, pconstr, ...){
	objfiter <- 0
	theta <- matrix(0, 7, dm + 2)
	dr1 <- pnew - par
	dr2 <- pnew - pold
	dr3 <- par - pold
	
	theta[1,2:(dm+2)] <- c(llp, par)
	theta[2,2:(dm+2)] <- c(llold, pold)
	theta[3,2:(dm+2)] <- c(llnew, pnew)
	temp <- search(pnew=pnew, dr=dr1, dm=dm, boundary=boundary, objfn=objfn, ...)
	objfiter <- objfiter + temp$objfiter
	theta[4,] <- temp$res
	dr1 <- temp$dr
	temp <- search(pnew=pnew, dr=dr2, dm=dm, boundary=boundary, objfn=objfn, ...)
	objfiter <- objfiter + temp$objfiter
	theta[5,] <- temp$res
	dr2 <- temp$dr
	temp <- search(pnew=pnew, dr=dr3, dm=dm, boundary=boundary, objfn=objfn, ...)
	theta[6,] <- temp$res
	objfiter <- objfiter + temp$objfiter
	
	if(prod(theta[4:6, 1])>0){
		a <- (theta[4,2]-llnew-theta[4,1]^2*(llp-llnew))/(theta[4,1]+theta[4,1]^2)
		c <- llp-llnew+a
		b <- (theta[5,2]-llnew-theta[5,1]^2*(llold-llnew))/(theta[5,1]+theta[5,1]^2)
		d <- llold-llnew+b
		e <- (theta[6,2]-llnew+theta[6,1]*(a-b)-theta[6,1]^2*(c+d))/(-2*theta[6,1]^2)
		#cat("a-e=", a,b,c,d,e,"\n")
		x <- -inv2(matrix(c(c,e,e,d),2,2))%*%c(a,b)/2
		#cat("x=", x,"\n")
		if(is.finite(prod(x))) {
			dr4 <- x[1]*dr1+x[2]*dr2
			temp <- search(pnew=pnew, dr=dr4, dm=dm, boundary=boundary, objfn=objfn, ...)
			theta[7,] <- temp$res
			objfiter <- objfiter + temp$objfiter
			sel <- which.max(theta[3:7,2])+2
		} else sel <- which.max(theta[3:6,2])+2
	} else sel <- which.max(theta[3:6,2])+2
	
	## Added if(pconstr...), JFB 13Feb2012
	if(pconstr(theta[sel,3:(dm+2)])) {
		pnew <- theta[sel,3:(dm+2)]
		llnew <- theta[sel,2]
	}
	#cat("llk=", theta[,2], "\n")
	#cat("sel=", sel, "\n")
	return(list(pnew=pnew, llnew=llnew, objfiter=objfiter))
}

########################################
########################################
## Main function for each algorithm ####
########################################
########################################
## input of function: should be a list containing the necessary components for running the algorithm
## output of function: should be a list containing "par", the current parameter value; "objfval", the current value of the objective function (=NULL if the algorithm does not evaluate the objective function)

bodyEM <- function(input.EM, ctrl, fixptfn, objfn, pconstr, project, iter, fpeval, objfeval, boundary, ...) {
	par <- fixptfn(input.EM$par, ...)
	fpeval <- fpeval + 1
	return(list(output=list(par=par, objfval=NULL), fpeval=fpeval, objfeval=objfeval))
}
bodyParaEM <- function(input.paraEM, ctrl, fixptfn, objfn, pconstr, project, iter, fpeval, objfeval, boundary, ...) {
	L2 <- -input.paraEM$objfval
	p0 <- input.paraEM$p0
	p1 <- input.paraEM$p1
	p2 <- input.paraEM$par
	
	if(ctrl$version=="geometric") {
		i <- 0
		t <- 1+ctrl$a^i*ctrl$h
	} else if(ctrl$version=="arithmetic") {
		i <- 1
		t <- 1 + i*ctrl$h
	} else { ## added 9Feb2012 (JFB)
		stop("for pem, version must be one of c('geometric','arithmetic')")	
	}
	pnew <- (1-t)^2*p0 + 2*t*(1-t)*p1 + t^2*p2
	no.pconstr <- !pconstr(pnew)
	if(no.pconstr)
		pnew <- project(pnew)
	Lnew <- try(-objfn(pnew, ...), silent=TRUE)
	objfeval <- objfeval + 1
	## Hui added !pconstr(pnew) here
	if(class(Lnew)=="try-error" | is.nan(Lnew) | no.pconstr | Lnew <= L2) { 
		p0 <- p2
		p1 <- fixptfn(p0, ...)
		p2 <- fixptfn(p1, ...)
		fpeval <- fpeval + 2
	} else {
		while(Lnew > L2) {
			pold <- pnew
			L2 <- Lnew
			i <- i+1
			t <- ifelse(ctrl$version=="arithmetic", 1 + i*ctrl$h, 1 + ctrl$a^i*ctrl$h)
			pnew <- (1-t)^2*p0 + 2*t*(1-t)*p1 + t^2*p2
			no.pconstr <- !pconstr(pnew)
			if(no.pconstr)
			   pnew <- project(pnew)
			Lnew <- try(-objfn(pnew, ...), silent=TRUE)
			objfeval <- objfeval + 1
			## Hui added !pconstr(pnew) [no.pconstr] here
			if(class(Lnew)=="try-error" | no.pconstr | is.nan(Lnew) | any(is.nan(unlist(pnew)))) {
				pnew <- pold
				Lnew <- L2
			}
		}
		p0 <- p1
		p1 <- p2
		p2 <- fixptfn(fixptfn(pold, ...), ...)
		fpeval <- fpeval + 2
	}
	
	L2 <- -objfn(p2, ...)
	objfeval <- objfeval + 1
	return(list(output=list(p0=p0, p1=p1, par=p2, objfval=-L2), fpeval=fpeval, objfeval=objfeval))
}
bodyDECME <- function(input.DECME, ctrl, fixptfn, objfn, pconstr, project, iter, fpeval, objfeval, boundary, ...) {
	pold <- input.DECME$pold
	par <- input.DECME$par
	llold <- -input.DECME$objfval.old
	llp <- -input.DECME$objfval
	
	pnew <- fixptfn(par, ...)
	fpeval <- fpeval + 1
	llnew <- try(-objfn(pnew, ...), silent=T)
	objfeval <- objfeval + 1
	
	if(ctrl$version %in% c(2,"v1")){ ## added "v1" 9Feb2012 (JFB)
		if(round((iter-1)/ctrl$dm)==(iter-1)/ctrl$dm) { ## restart
			res <- line.search(pnew, pnew, par, llnew, ctrl$tol_op, boundary, objfn, ...) 	# restart
			objfeval <- objfeval + res$optimiter
		} else {
			res <- line.search(pnew, pnew, par, llnew, ctrl$tol_op, boundary, objfn, ...)
			objfeval <- objfeval + res$optimiter
			res <- line.search(res$pnew, res$pnew, pold, res$llnew, ctrl$tol_op, boundary, objfn, ...)
			objfeval <- objfeval + res$optimiter
		}
		pnew <- res$pnew
		llnew <- res$llnew
	} else if(ctrl$version=="2s") {
		if(iter > 1) {
			res <- optimize.2d.app(pnew, pold, par, llnew, llold, llp, dm=ctrl$dm, boundary, objfn, pconstr, ...)
			objfeval <- objfeval + res$objfiter
			pnew <- res$pnew
			llnew <- res$llnew
		}
	} else if(ctrl$version=="v2") {
		if(iter > 1) {
			res <- line.search(pnew, pnew, pold, llnew, ctrl$tol_op, boundary, objfn, ...) 	# restart
			objfeval <- objfeval + res$optimiter
			pnew <- res$pnew
			llnew <- res$llnew
		}
	} else if(ctrl$version=="v3") {
		if(iter > 1) {
			res <- line.search(pnew, par, pold, llnew, ctrl$tol_op, boundary, objfn, ...) 	# restart
			objfeval <- objfeval + res$optimiter
			pnew <- res$pnew
			llnew <- res$llnew
		}
	} else { ## added 9Feb2012 (JFB)
		stop("For decme, version must be one of c('2','2s','v2','v3')")	
	}
	
	return(list(output=list(pold=par, par=pnew, objfval.old=-llp, objfval=-llnew), fpeval=fpeval, objfeval=objfeval))
}
bodyQuasiNewton <- function(input.quasiNewton, ctrl, fixptfn, objfn, pconstr, project, iter, fpeval, objfeval, boundary, ...) {
	par   <- input.quasiNewton$par
	pnew  <- input.quasiNewton$pnew
	pnew2 <- input.quasiNewton$pnew2
	U     <- input.quasiNewton$U
	V     <- input.quasiNewton$V
	
	
	pold <- par
	##tmp <- try(solve(t(U)%*%U-t(U)%*%V),silent=TRUE)
	tmp <- try(solve(crossprod(U, U-V)),silent=TRUE)
	if(class(tmp)=='try-error')
		# reassign p.next to be
		p.next <- pnew2
	else
	{
		tmp1 <- V %*% tmp
		tmp2 <- crossprod(U, par - pnew)
		p.next <- as.vector( pnew - tmp1 %*% tmp2 )
		if(!pconstr(pnew))
			pnew <- project(pnew)
	}
	
	if(!pconstr(p.next))
		p.next <- pnew2
	if(!is.null(objfn))
	{
		lnew <- objfn(p.next, ...)
		objfeval <- objfeval + 1
		if(any(p.next!=pnew2))
		{
			l2 <- objfn(pnew2, ...)
			objfeval <- objfeval + 1
			if(l2<lnew)
			{
				p.next <- pnew2
				lnew   <- l2
			}
		}
	}
	
	par <- p.next
	pnew <- fixptfn(par, ...)
	fpeval <- fpeval + 1
	U <- cbind(U[,-1],pnew-par)
	pnew2 <- fixptfn(pnew, ...)
	pfeval <- fpeval + 1
	V <- cbind(V[,-1],pnew2-pnew)
	if(!is.null(objfn))
		return(list(output=list(par=par, pnew=pnew, pnew2=pnew2, V=V, U=U, objfval=lnew), fpeval=fpeval, objfeval=objfeval))
	
	return(list(output=list(par=par, pnew=pnew, pnew2=pnew2, V=V, U=U, objfval=NULL), fpeval=fpeval, objfeval=objfeval))
}
bodySquarem1 <- function(input.squarem1, ctrl, fixptfn, objfn, pconstr, project, iter, fpeval, objfeval, boundary, ...) {
	p <- input.squarem1$par
	lold <- input.squarem1$objfval
	step.min <- input.squarem1$step.min
	step.max <- input.squarem1$step.max
	
	extrap <- TRUE
	p1 <- try(fixptfn(p, ...), silent = TRUE)
	fpeval <- fpeval + 1
	if (class(p1) == "try-error" | any(is.nan(unlist(p1)))) 
		stop("Error in function evaluation")
	q1 <- p1 - p
	sr2 <- crossprod(q1)
	p2 <- try(fixptfn(p1, ...), silent = TRUE)
	pfeval <- fpeval + 1
	if (class(p2) == "try-error" | any(is.nan(unlist(p2)))) 
		stop("Error in function evaluation")
	q2 <- p2 - p1
	sq2 <- sqrt(crossprod(q2))
	sv2 <- crossprod(q2 - q1)
	srv <- crossprod(q1, q2 - q1)
	alpha <- switch(ctrl$version, -srv/sv2, -sr2/srv, sqrt(sr2/sv2))
	alpha <- max(step.min, min(step.max, alpha))
	p.new <- p + 2 * alpha * q1 + alpha^2 * (q2 - q1)
	if(!pconstr(p.new))
		p.new <- project(p.new)
	# Modified (JFB 30Jan2012)
	# if (abs(alpha - 1) > 0.01) {
	if (isTRUE(abs(alpha - 1) > 0.01) ) {
		p.new <- try(fixptfn(p.new, ...), silent = TRUE)
		fpeval <- fpeval + 1
	}
	if (class(p.new) == "try-error" | any(is.nan(p.new)) | !pconstr(p.new)) {
		p.new <- p2
		lnew <- try(objfn(p2, ...), silent = TRUE)
		objfeval <- objfeval + 1
		# Modified (JFB 30Jan2012)
		# if (alpha == step.max) 
		if (isTRUE(all.equal(alpha, step.max)))
			step.max <- max(ctrl$step.max0, step.max/ctrl$mstep)
		alpha <- 1
		extrap <- FALSE
	} else {
		if (is.finite(ctrl$objfn.inc)) {
			lnew <- try(objfn(p.new, ...), silent = TRUE)
			objfeval <- objfeval + 1
		}
		else lnew <- lold
		if (class(lnew) == "try-error" | is.nan(lnew) | (lnew > 
					lold + ctrl$objfn.inc)) {
			p.new <- p2
			lnew <- try(objfn(p2, ...), silent = TRUE)
			objfeval <- objfeval + 1
			# Modified (JFB 30Jan2012)
			# if (alpha == step.max) 
			if (isTRUE(all.equal(alpha, step.max)))
				step.max <- max(ctrl$step.max0, step.max/ctrl$mstep)
			alpha <- 1
			extrap <- FALSE
		}
	}
	# Modified (JFB 30Jan2012)
	# if (alpha == step.max) 
	if (isTRUE(all.equal(alpha, step.max)))
		step.max <- ctrl$mstep * step.max
	# Modified (JFB 30Jan2012)
	# if (step.min < 0 & alpha == step.min) 
	if (step.min < 0 & isTRUE(all.equal(alpha, step.min)))
		step.min <- ctrl$mstep * step.min
	p <- p.new
	if (!is.nan(lnew)) 
		lold <- lnew
	
	return(list(output=list(par=p, objfval=lold, step.min=step.min, step.max=step.max), fpeval=fpeval, objfeval=objfeval))
}
bodySquarem2 <- function(input.squarem2, ctrl, fixptfn, objfn, pconstr, project, iter, fpeval, objfeval, boundary, ...) {
	p <- input.squarem2$par
	step.min <- input.squarem2$step.min
	step.max <- input.squarem2$step.max
	
	extrap <- TRUE
	p1 <- try(fixptfn(par, ...), silent = TRUE)
	fpeval <- fpeval + 1
	if (class(p1) == "try-error" | any(is.nan(unlist(p1)))) 
		break
	q1 <- p1 - par
	sr2 <- crossprod(q1)
	# if (sqrt(sr2) < tol) 
	# break
	p2 <- try(fixptfn(p1, ...), silent = TRUE)
	fpeval <- fpeval + 1
	if (class(p2) == "try-error" | any(is.nan(unlist(p2)))) 
		break
	q2 <- p2 - p1
	sq2 <- sqrt(crossprod(q2))
	res <- sq2
	# if (sq2 < tol) 
	# break
	sv2 <- crossprod(q2 - q1)
	srv <- crossprod(q1, q2 - q1)
	alpha <- switch(ctrl$version, -srv/sv2, -sr2/srv, sqrt(sr2/sv2))
	alpha <- max(step.min, min(step.max, alpha))
	p.new <- par + 2 * alpha * q1 + alpha^2 * (q2 - q1)
	if(!pconstr(p.new))
		p.new <- project(p.new)
	# Modified (JFB 30Jan2012)
	# if (abs(alpha - 1) > 0.01) {
	if (isTRUE(abs(alpha - 1) > 0.01) ) {	
		ptmp <- try(fixptfn(p.new, ...), silent = TRUE)
		fpeval <- fpeval + 1
		if (class(ptmp) == "try-error" | any(is.nan(unlist(ptmp))) | !pconstr(ptmp)) {
			p.new <- p2
			# Modified (JFB 30Jan2012)
			# if (alpha == step.max) 
			if (isTRUE(all.equal(alpha, step.max)))
				step.max <- max(ctrl$step.max0, step.max/ctrl$mstep)
			alpha <- 1
			extrap <- FALSE
		}
		else {
			res <- sqrt(crossprod(ptmp - p.new))
			parnorm <- sqrt(crossprod(p2)/length(p2))
			kres <- ctrl$kr * (1 + parnorm) + sq2
			p.new <- if (res <= kres) 
						ptmp
					else p2
			if (res > kres) {
				# Modified (JFB 30Jan2012)
				# if (alpha == step.max) 
				if (isTRUE(all.equal(alpha, step.max)))
					step.max <- max(ctrl$step.max0, step.max/ctrl$mstep)
				alpha <- 1
				extrap <- FALSE
			}
		}
	}
	# Modified (JFB 30Jan2012)
	# if (alpha == step.max) 
	if (isTRUE(all.equal(alpha, step.max)))
		step.max <- ctrl$mstep * step.max
	# Modified (JFB 30Jan2012)
	# if (step.min < 0 & alpha == step.min) 
	if (step.min < 0 & alpha == step.min)
		step.min <- ctrl$mstep * step.min
	par <- p.new
	extrap <- TRUE
	p1 <- try(fixptfn(par, ...), silent = TRUE)
	fpeval <- fpeval + 1
	if (class(p1) == "try-error" | any(is.nan(unlist(p1)))) 
		break
	q1 <- p1 - par
	sr2 <- crossprod(q1)
	# if (sqrt(sr2) < tol) 
	# break
	p2 <- try(fixptfn(p1, ...), silent = TRUE)
	fpeval <- fpeval + 1
	if (class(p2) == "try-error" | any(is.nan(unlist(p2)))) 
		break
	q2 <- p2 - p1
	sq2 <- sqrt(crossprod(q2))
	res <- sq2
	# if (sq2 < tol) 
	# break
	sv2 <- crossprod(q2 - q1)
	srv <- crossprod(q1, q2 - q1)
	alpha <- switch(ctrl$version, -srv/sv2, -sr2/srv, sqrt(sr2/sv2))
	alpha <- max(step.min, min(step.max, alpha))
	p.new <- par + 2 * alpha * q1 + alpha^2 * (q2 - q1)
	if(!pconstr(p.new))
		p.new <- project(p.new)
	# Modified (JFB 30Jan2012)
	# if (abs(alpha - 1) > 0.01) {
	if (isTRUE(abs(alpha - 1) > 0.01) ) {
		ptmp <- try(fixptfn(p.new, ...), silent = TRUE)
		fpeval <- fpeval + 1
		if (class(ptmp) == "try-error" | any(is.nan(unlist(ptmp))) | !pconstr(ptmp)) {
			p.new <- p2
			if (alpha == step.max) 
				step.max <- max(ctrl$step.max0, step.max/ctrl$mstep)
			alpha <- 1
			extrap <- FALSE
		}
		else {
			res <- sqrt(crossprod(ptmp - p.new))
			parnorm <- sqrt(crossprod(p2)/length(p2))
			kres <- ctrl$kr * (1 + parnorm) + sq2
			p.new <- if (res <= kres) 
						ptmp
					else p2
			if (res > kres) {
				if (alpha == step.max) 
					step.max <- max(ctrl$step.max0, step.max/ctrl$mstep)
				alpha <- 1
				extrap <- FALSE
			}
		}
	}
	if (alpha == step.max) 
		step.max <- ctrl$mstep * step.max
	if (step.min < 0 & alpha == step.min) 
		step.min <- ctrl$mstep * step.min
	par <- p.new
	
	return(list(output=list(par=par, objfval=NULL, step.min=step.min, step.max=step.max), fpeval=fpeval, objfeval=objfeval))
}
bodyCyclem1 <- function(input.cyclem1, ctrl, fixptfn, objfn, pconstr, project, iter, fpeval, objfeval, boundary, ...) {
	par <- input.cyclem1$par
	lold <- input.cyclem1$objfval
	
	K <- ctrl$K
	
	extrap <- TRUE
	p <- matrix(NA, nrow = length(par), ncol = K + 2)
	U <- matrix(NA, nrow = length(par), ncol = K + 1)
	p[, 1] <- par
	for (i in 2:(K + 2)) {
		p[, i] <- try(fixptfn(p[, i - 1], ...), silent = TRUE)
		fpeval <- fpeval + 1
		if (class(p[, i]) == "try-error" | any(is.nan(unlist(p[, 
										i]))) | any(!is.numeric(p[, i]))) 
			stop("Error in function evaluation \n")
		U[, i - 1] <- p[, i] - p[, i - 1]
	}
	sqKp1 <- sqrt(c(crossprod(U[, K + 1])))
	if (ctrl$version == "rre") {
		coef <- try(solve(qr(crossprod(U), LAPACK = TRUE, 
								tol = 1e-14), rep(1, K + 1)), silent = TRUE)
		if (class(coef) == "try-error" | any(is.nan(coef))) {
			extrap <- FALSE
		}
		else {
			# Modified (JFB 30Jan2012)
			# if (abs(sum(coef)) < 1e-07) 
			if (isTRUE(all.equal(abs(sum(coef)), 1.e-07, tolerance =  1.e-07)))
				extrap <- FALSE
			coef <- coef/sum(coef)
		}
	}
	if (ctrl$version == "mpe") {
		coef <- try(solve(qr(U[, -(K + 1)], LAPACK = TRUE, 
								tol = 1e-14), -U[, K + 1]), silent = TRUE)
		if (class(coef) == "try-error" | any(is.nan(coef))) {
			extrap <- FALSE
		}
		else {
			coef <- c(coef, 1)
			# Modified (JFB 30Jan2012)
			# if (abs(sum(coef)) < 1e-07) 
			if (isTRUE(all.equal(abs(sum(coef)), 1.e-07, tolerance =  1.e-07)))
				extrap <- FALSE
			coef <- coef/sum(coef)
		}
	}
	if (!extrap) {
		par <- p[, ncol(p)]
		res <- sqKp1
		##next
		return(list(output=list(par=par, objfval=lold), fpeval=fpeval, objfeval=objfeval))
	}
	if (!ctrl$square) 
	{
		pnew <- c(p[, -(K + 2)] %*% coef)
		if(!pconstr(pnew))
			pnew <- project(pnew)
	}
	if (ctrl$square) {
		pnew <- rep(0, length(par))
		if (K > 1) {
			for (i in 1:(K - 1)) {
				p <- try(cbind(p, fixptfn(p[, i + K + 1], ...)), silent = TRUE)
				fpeval <- fpeval + 1
				if (class(p) == "try-error" | any(is.nan(unlist(p))) | any(!is.numeric(p))) 
					stop("Error in function evaluation \n")
			}
		}
		for (i in 0:K) {
			for (j in 0:K) {
				pnew <- pnew + coef[i + 1] * coef[j + 1] * p[, i + j + 1]
			}
		}
		if(!pconstr(pnew))
			pnew <- project(pnew)
		sqKp1 <- sqrt(c(crossprod(p[, 2 * K] - p[, 2 * K - 
												1])))
	}
	if (extrap) {
		ptmp <- try(fixptfn(pnew, ...), silent = TRUE)
		fpeval <- fpeval + 1
		res <- sqrt(c(crossprod(ptmp - pnew)))
	}
	if (class(ptmp) == "try-error" | any(is.nan(unlist(ptmp))) | !pconstr(ptmp)) {
		pnew <- p[, ncol(p)]
		res <- sqKp1
		extrap <- FALSE
	}
	else pnew <- ptmp
	lnew <- try(objfn(pnew, ...), silent = TRUE)
	objfeval <- objfeval + 1
	if (class(lnew) == "try-error" | is.nan(lnew) | (lnew > 
				lold + ctrl$objfn.inc)) {
		pnew <- p[, ncol(p)]
		res <- sqKp1
		extrap <- FALSE
	}
	else lold <- lnew
	par <- pnew
	
	return(list(output=list(par=par, objfval=lold), fpeval=fpeval, objfeval=objfeval))
}
bodyCyclem2 <- function(input.cyclem2, ctrl, fixptfn, objfn, pconstr, project, iter, fpeval, objfeval, boundary, ...) {
	par <- input.cyclem2$par
	
	K <- ctrl$K
	
	extrap <- TRUE
	p <- matrix(NA, nrow = length(par), ncol = K + 2)
	U <- matrix(NA, nrow = length(par), ncol = K + 1)
	p[, 1] <- par
	for (i in 2:(K + 2)) {
		p[, i] <- try(fixptfn(p[, i - 1], ...), silent = TRUE)
		fpeval <- fpeval + 1
		if (class(p[, i]) == "try-error" | any(is.nan(unlist(p[, 
										i]))) | any(!is.numeric(p[, i]))) 
			stop("Error in function evaluation")
		U[, i - 1] <- p[, i] - p[, i - 1]
	}
	sqKp1 <- sqrt(c(crossprod(U[, K + 1])))
	if (ctrl$version == "rre") {
		coef <- try(solve(qr(crossprod(U), LAPACK = TRUE, 
								tol = 1e-14), rep(1, K + 1)), silent = TRUE)
		if (class(coef) == "try-error" | any(is.nan(coef))) {
			extrap <- FALSE
		}
		else {
			# Modified (JFB 30Jan2012)
			# if (abs(sum(coef)) < 1e-07) 
			if (isTRUE(all.equal(abs(sum(coef)), 1.e-07, tolerance = 1.e-07)))
				extrap <- FALSE
			coef <- coef/sum(coef)
		}
	}
	if (ctrl$version == "mpe") {
		coef <- try(solve(qr(U[, -(K + 1)], LAPACK = TRUE, 
								tol = 1e-14), -U[, K + 1]), silent = TRUE)
		if (class(coef) == "try-error" | any(is.nan(coef))) {
			extrap <- FALSE
		}
		else {
			coef <- c(coef, 1)
			# Modified (JFB 30Jan2012)
			# if (abs(sum(coef)) < 1e-07) 
			if (isTRUE(all.equal(abs(sum(coef)), 1.e-07, tolerance = 1.e-07)))
				extrap <- FALSE
			coef <- coef/sum(coef)
		}
	}
	if (!extrap) {
		par <- p[, ncol(p)]
		res <- sqKp1
		return(list(output=list(par=par, objfval=NULL), fpeval=fpeval, objfeval=objfeval))
	}
	if (!ctrl$square) 
	{
		pnew <- c(p[, -(K + 2)] %*% coef)
		if(!pconstr(pnew))
			pnew <- project(pnew)
	}
	if (ctrl$square) {
		pnew <- rep(0, length(par))
		if (K > 1) {
			for (i in 1:(K - 1)) {
				p <- try(cbind(p, fixptfn(p[, i + K + 1], ...)), 
						silent = TRUE)
				fpeval <- fpeval + 1
				if (class(p) == "try-error" | any(is.nan(unlist(p))) | 
						any(!is.numeric(p))) 
					stop("Error in function evaluation")
			}
		}
		for (i in 0:K) {
			for (j in 0:K) {
				pnew <- pnew + coef[i + 1] * coef[j + 1] * 
						p[, i + j + 1]
			}
		}
		if(!pconstr(pnew))
			pnew <- project(pnew)
		sqKp1 <- sqrt(c(crossprod(p[, 2 * K] - p[, 2 * K - 
												1])))
	}
	ptemp <- try(fixptfn(pnew, ...), silent = TRUE)
	fpeval <- fpeval + 1
	if (class(ptemp) == "try-error" | any(is.nan(ptemp)) | !pconstr(ptemp)) {
		par <- p[, ncol(p)]
		res <- sqKp1
		extrap <- FALSE
	}
	else {
		res <- sqrt(c(crossprod(ptemp - pnew)))
		parnorm <- as.numeric(sqrt(crossprod(p[, ncol(p)])/length(par)))
		if (res <= (ctrl$kr * (1 + parnorm) + sqKp1)) {
			par <- ptemp
			extrap <- TRUE
		}
		else {
			par <- p[, ncol(p)]
			res <- sqKp1
			extrap <- FALSE
		}
	}
	
	return(list(output=list(par=par, objfval=NULL), fpeval=fpeval, objfeval=objfeval))
}

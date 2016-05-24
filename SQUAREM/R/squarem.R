squarem <- function(par, fixptfn, objfn, ... , control=list()) {
#
# The Master function for SQUAREM acceleration of "ANY" fixed-point iteration
# A wrapper that calls appropriate function below (squarem1, squarem2, cyclem1, cyclem2)
#
# Author:  Ravi Varadhan, Johns Hopkins University
#
# See Varadhan and Roland (Scandinavian J Statistics 2008)
# See Roland, Varadhan and Frangakis (Applied Numerical Mathematics 2007)
#
# Last modified: June 25, 2010
# Last modified: August 09, 2010
# Last modified: December 31, 2010
#######################################################################
#
# par = starting value of parameter vector
# fixptfn = fixed-point iteration F(x)
# for which the solution: F(x*) = x* is sought
# objfn = underlying objective function which is minimized at x*
# method = 1, 2, or 3, indicating the type of steplength to be used
#
#
control.default <- list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,
	tol=1.e-07, maxiter=1500, trace=FALSE)
namc <- names(control)
if (!all(namc %in% names(control.default))) 
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])

ctrl <- modifyList(control.default, control)

if (ctrl$K > 1 & !(ctrl$method %in% c("rre", "mpe")) ) ctrl$method <- "rre"
if (ctrl$K == 1 & !(ctrl$method %in% c(1,2,3)) ) ctrl$method <- 3

if (!missing(objfn) ) {
if (ctrl$K == 1) sqobj <- squarem1(par, fixptfn, objfn, ... , control=ctrl)
else if (ctrl$K > 1 | ctrl$method %in% c("rre", "mpe")) sqobj <- cyclem1(par, fixptfn, objfn, ... , control=ctrl)
} else {
if (ctrl$K == 1) sqobj <- squarem2(par, fixptfn, ... , control=ctrl)
else if (ctrl$K > 1 | ctrl$method %in% c("rre", "mpe"))  sqobj <- cyclem2(par, fixptfn, ... , control=ctrl)
}
return(sqobj)
}

#######################################################################
# Partially-monotone, globally-convergent SQUAREM acceleration
# Can accelerate EM, MM, and other fixed-point iterations
#  Here we use an adaptive maximum step-size strategy
#  This strategy is effective and represents a good trade-off between speed and stability
# See Varadhan and Roland (Scandinavian J Statistics 2008)
#
# Date: August 22, 2007
# Adaptive step-max strategy modified:  May 2008
# Last modified:  November 2008
# Last modified: October 2009
# Last modified: Feb 23 2010
#######################################################################
squarem1 <- function(par, fixptfn, objfn, ... , control=list()) {
# par = starting value of parameter vector
# fixptfn = fixed-point iteration F(x)
# for which the solution: F(x*) = x* is sought
# objfn = underlying objective function which is minimized at x*
#
#
control.default <- list(K=1, square=TRUE, method=3, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1, 
	tol=1.e-07, maxiter=1500, trace=FALSE)

namc <- names(control)
if (!all(namc %in% names(control.default))) 
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
ctrl <- modifyList(control.default, control)

#####
# method = 1, 2, or 3, indicating the type of steplength to be used
# A key parameter is `step.min0'. This can be either positive or negative if the eigenvalues of the Jacobian of fixed-point mapping are all positive;
# We set it to +1 for contraction mappings such as EM and MM algorithms
# Must be negative, e.g. -1, if the fixed-point mapping can have negative eiegnvalues
#####
# parameter "objfn.inc" dictates the amount of non-monotonicity in the objective function 
# setting objfn.inc=0 would enforce monotonicity, whereas objfn.inc=Inf would be a non-monotonic scheme
# The defalut objfn.inc=1 would enforce monotonicity far from solution, but allows for non-monotonicity closer to solution
#
    method <- ctrl$method
    maxiter <- ctrl$maxiter
    tol <- ctrl$tol
    step.min <- ctrl$step.min0
    step.max0 <- ctrl$step.max0
    step.max <- ctrl$step.max0
    mstep <- ctrl$mstep
    objfn.inc <- ctrl$objfn.inc
    trace <- ctrl$trace

if (trace) cat("Squarem-1 \n")

if (missing(objfn)) stop("\n squarem2 should be used if objective function is not available \n\n")

iter <- 1
objval <- rep(NA,1)
p <- par

lold <- objfn(p, ...)
leval <- 1
if (trace) cat("Objective fn: ", lold, "\n")
feval <- 0
conv <- TRUE

while (feval < maxiter) {
	extrap <- TRUE
	p1 <- try(fixptfn(p, ...),silent=TRUE)
	feval <- feval + 1
	if (class(p1) == "try-error" | any(is.nan(unlist(p1)))) stop("Error in function evaluation")
	q1 <- p1 - p
	sr2 <- crossprod(q1)
	if (sqrt(sr2) < tol) break

	p2 <- try(fixptfn(p1, ...),silent=TRUE)
	feval <- feval + 1
	if (class(p2) == "try-error" | any(is.nan(unlist(p2)))) stop("Error in function evaluation")

	q2 <- p2 - p1
	sq2 <- sqrt(crossprod(q2))
	if (sq2 < tol) break
	sv2 <- crossprod(q2-q1)
	srv <- crossprod(q1, q2-q1)

	alpha <- switch(method, -srv/sv2, -sr2/srv, sqrt(sr2/sv2)) 

 	alpha <- max(step.min, min(step.max, alpha))
	p.new <- p + 2*alpha*q1 + alpha^2*(q2-q1)
	if (abs(alpha - 1) > 0.01 ) {
		p.new <- try(fixptfn(p.new , ...),silent=TRUE)   # stabilization step
		feval <- feval + 1
	}

	if (class(p.new) == "try-error" | any(is.nan(p.new))) {
	p.new <- p2
	lnew <- try(objfn(p2, ...), silent=TRUE)
	leval <- leval + 1
	if (alpha == step.max) step.max <- max(step.max0, step.max/mstep)
	alpha <- 1
	extrap <- FALSE
	} else {
		if (is.finite(objfn.inc)) {
			lnew <- try(objfn(p.new, ...), silent=TRUE)
			leval <- leval + 1
		} else lnew <- lold
		if (class(lnew) == "try-error" | is.nan(lnew) | 
		(lnew > lold + objfn.inc)) {
			p.new <- p2
			lnew <- try(objfn(p2, ...), silent=TRUE)
			leval <- leval + 1
			if (alpha==step.max) step.max <- max(step.max0, step.max/mstep)
			alpha <- 1
			extrap <- FALSE
		}
	}	

	if (alpha == step.max) step.max <- mstep*step.max
	if (step.min < 0 & alpha == step.min) step.min <- mstep*step.min

	p <- p.new
	if (!is.nan(lnew)) lold <- lnew 
	if (trace) cat("Objective fn: ", lnew, "  Extrapolation: ", extrap, "  Steplength: ", alpha, "\n")
	iter <- iter+1
	
}
	if (feval >= maxiter) conv <- FALSE
	if (is.infinite(objfn.inc)) {
		lold <- objfn(p, ...)
		leval <- leval + 1
		}


return(list(par=p, value.objfn=lold, iter= iter, fpevals=feval, objfevals = leval, convergence=conv))

}
###################################################################
squarem2 <- function(par, fixptfn, ... , control=list() ) {

# par = starting value of parameter vector
# fixptfn = fixed-point iteration F(x)
# for which the solution: F(x*) = x* is sought
# method = 1, 2, or 3, indicating the type of steplength to be used
#
# See Varadhan and Roland (Scandinavian J Statistics 2008)
#
control.default <- list(K=1, square=TRUE, method=3, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1, 
	tol=1.e-07, maxiter=1500, trace=FALSE)

namc <- names(control)
if (!all(namc %in% names(control.default))) 
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
ctrl <- modifyList(control.default, control)

#####
# method = 1, 2, or 3, indicating the type of steplength to be used
# A key parameter is `step.min0'. This can be either positive or negative if the eigenvalues of the Jacobian of fixed-point mapping are all positive;
# We set it to +1 for contraction mappings such as EM and MM algorithms
# Must be negative, e.g. -1, if the fixed-point mapping can have negative eiegnvalues
#####
# Parameter "kr" dictates the amount of non-monotonicity in the fixed-point residual 
# setting kr=0 would try to enforce monotonicity; kres=Inf would be a non-monotonic scheme
#
    method <- ctrl$method
    maxiter <- ctrl$maxiter
    tol <- ctrl$tol
    kr <- ctrl$kr
    step.min <- ctrl$step.min0
    step.max0 <- ctrl$step.max0
    step.max <- ctrl$step.max0
    mstep <- ctrl$mstep
    trace <- ctrl$trace
if (trace) cat("Squarem-2 \n")

iter <- 1
feval <- 0
kount <- 0
conv <- TRUE

while (feval < maxiter) {
	extrap <- TRUE
	p1 <- try(fixptfn(par, ...),silent=TRUE)
	feval <- feval + 1
	if (class(p1) == "try-error" | any(is.nan(unlist(p1)))) break
	q1 <- p1 - par
	sr2 <- crossprod(q1)
	if (sqrt(sr2) < tol) break

	p2 <- try(fixptfn(p1, ...),silent=TRUE)
	feval <- feval + 1
	if (class(p2) == "try-error" | any(is.nan(unlist(p2)))) break
	q2 <- p2 - p1
	sq2 <- sqrt(crossprod(q2))
	res <- sq2
	if (sq2 < tol) break
	sv2 <- crossprod(q2-q1)
	srv <- crossprod(q1, q2-q1)	

	alpha <- switch(method, -srv/sv2, -sr2/srv, sqrt(sr2/sv2)) 

	alpha <- max(step.min, min(step.max, alpha))
	p.new <- par + 2*alpha*q1 + alpha^2*(q2-q1) 

	if (abs(alpha - 1) > 0.01 ) {
	ptmp <- try(fixptfn(p.new, ...),silent=TRUE)   # stabilization step
	feval <- feval + 1
		if (class(ptmp) == "try-error" | any(is.nan(unlist(ptmp))) ) {
		p.new <- p2
		if (alpha == step.max) step.max <- max(step.max0, step.max/mstep)
		alpha <- 1 
		extrap <- FALSE
      		}  else {
			res <- sqrt(crossprod(ptmp - p.new))
			parnorm <- sqrt(crossprod(p2)/length(p2))
			kres <- kr * (1 + parnorm) + sq2 
			p.new <- if (res <= kres) ptmp else p2
			if (res > kres) {
				if (alpha == step.max) step.max <- max(step.max0, step.max/mstep)
				alpha <- 1
				extrap <- FALSE
			}
		}
	}

	if (alpha == step.max) step.max <- mstep*step.max
	if (step.min < 0 & alpha == step.min) step.min <- mstep*step.min

	if (trace) cat("Residual: ", res, "  Extrapolation: ", extrap, "  Steplength: ", alpha, "\n")

	par <- p.new
	iter <- iter+1
}
	if (feval >= maxiter) conv <- FALSE

return(list(par=par, value.objfn=NA, iter = iter, fpevals=feval, objfevals = 0, convergence=conv))

}
#######################################################################
# Squared-Cycled vector extrapolation for accelerating vector sequences
# See Roland, Varadhan and Frangakis (Applie Numerical Mathematics 2007)
# Author:  Ravi Varadhan, Johns Hopkins University
# Latest modofications:: Feb 23, 2010
# 
###########################
cyclem1 <- function(par, fixptfn, objfn, control=list(), ...) {
#
# par = starting value of parameter vector
# fixptfn = fixed-point iteration F(x)
# for which the solution: F(x*) = x* is sought
# objfn = underlying objective function which is minimized at x*
# 
control.default <- list(K=2, method="rre", square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1, 
	tol=1.e-07, maxiter=1500, trace=FALSE)

namc <- names(control)
if (!all(namc %in% names(control.default))) 
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
ctrl <- modifyList(control.default, control)

#
# method = reduced-rank or minimal-polynomial extrapolation
# K = order of extrapolation scheme; K=2,3,4 are typical choices
# square = a logical variable indicating whether or not "squaring" is used

    method <- ctrl$method
    K <- ctrl$K
    square <- ctrl$square
    tol <- ctrl$tol
    maxiter <- ctrl$maxiter
    objfn.inc <- ctrl$objfn.inc
    trace <- ctrl$trace
 if (trace) cat("Cyclem-1 \n")

if (missing(objfn)) stop("\n cyclem2 should be used if objective function is not available \n\n")

if ( K > length(par)) {
cat("K is too large.  Decrease it.  \n")
return()
}

iter <- 0
lold <- objfn(par, ...)
leval <- 1

ptmp <- fixptfn(par,...)
res <- sqrt(c(crossprod(ptmp - par)))
par <- ptmp
feval <- 1
conv <- TRUE

while (feval < maxiter & res > tol) {
	iter <- iter+1
	extrap <- TRUE
	p <- matrix(NA, nrow=length(par), ncol=K+2)
	U <- matrix(NA, nrow=length(par), ncol=K+1)
	p[,1] <- par
	for (i in 2:(K+2)){
	p[,i] <- try(fixptfn(p[,i-1], ...), silent=TRUE)
	if (class(p[,i]) == "try-error" | any(is.nan(unlist(p[,i]))) | any(!is.numeric(p[,i]))) stop("Error in function evaluation \n")
	U[,i-1] <- p[, i] - p[, i-1]
	}
	feval <- feval + K + 1 
	sqKp1 <- sqrt(c(crossprod(U[, K+1])))

	if (method == "rre"){
	coef <- try(solve(qr(crossprod(U), LAPACK=TRUE, tol=1.e-14), rep(1,K+1)), silent=TRUE)
	if (class(coef) == "try-error" | any(is.nan(coef)) ) {
		extrap <- FALSE
	} else {
		if (abs(sum(coef)) < 1.e-07) extrap <- FALSE
		coef <- coef/sum(coef)
		}
	}
	if (method == "mpe"){
	coef <- try(solve(qr(U[,-(K+1)], LAPACK=TRUE, tol=1.e-14), -U[,K+1]), silent=TRUE)
	if (class(coef) == "try-error" | any(is.nan(coef))) {
		extrap <- FALSE
	} else {
		coef <- c(coef, 1)
		if (abs(sum(coef)) < 1.e-07) extrap <- FALSE
		coef <- coef/sum(coef)
		}
	}

	if (!extrap) {
		par <- p[, ncol(p)]
		res <- sqKp1
		iter <- iter + 1
 		if (trace) cat(" objective fn: ", lold, "extrapolation: ", extrap, "\n")
		next 
	}

	if (!square) pnew <- c(p[,-(K+2)] %*% coef)
	if (square) {
		pnew <- rep(0, length(par))
		if (K > 1) {
		for (i in 1:(K-1)) {
			p <- try(cbind(p, fixptfn(p[,i+K+1],...)), silent=TRUE)
			if (class(p) == "try-error" | any(is.nan(unlist(p))) | any(!is.numeric(p)))  stop("Error in function evaluation \n")
			}
		feval <- feval + K - 1
		}
		for (i in 0:K){
		for (j in 0:K){
			pnew <- pnew + coef[i+1]* coef[j+1] * p[,i+j+1]
		} }
		sqKp1 <- sqrt(c(crossprod(p[, 2*K] - p[, 2*K-1])))
	}

#  Perform one stabilization step between cycles
		if (extrap) {
		ptmp <- try(fixptfn(pnew, ...), silent=TRUE)
		res <- sqrt(c(crossprod(ptmp - pnew)))
		feval <- feval + 1
		} 
	if (class(ptmp) == "try-error" | any(is.nan(unlist(ptmp))) ) {
		pnew <- p[, ncol(p)]
		res <- sqKp1
		extrap <- FALSE
		} else pnew <- ptmp
		
	lnew <- try(objfn(pnew, ...), silent=TRUE)
	leval <- leval + 1
	if (class(lnew) == "try-error" | is.nan(lnew) | (lnew > lold + objfn.inc)) {
		pnew <- p[, ncol(p)]
		res <- sqKp1
		extrap <- FALSE
 		if (trace) cat(" objective fn: ", lold, "extrapolation: ", extrap, "\n")
	} else lold <- lnew

	par <- pnew
 	if (trace) cat(" objective fn: ", lold, "extrapolation: ", extrap, "\n")
}  # Main loop complete

	if (feval >= maxiter) conv <- FALSE

return(list(par=par, value.objfn=lold, iter = iter, fpevals=feval, objfevals = leval, convergence=conv))

}
###################################################################
# Squared-Cycled vector extrapolation for accelerating vector sequences
# See Roland, Varadhan and Frangakis (Applie Numerical Mathematics 2007)
# Author:  Ravi Varadhan, Johns Hopkins University
# Latest modofications:: Feb 23, 2010
# 
###########################
cyclem2 <- function(par, fixptfn, control=list(), ...) {
# par = starting value of parameter vector
# fixptfn = fixed-point iteration F(x)
# for which the solution: F(x*) = x* is sought
# method = reduced-rank or minimal-polynomial extrapolation
#
control.default <- list(K=2, method="rre", square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1, 
	tol=1.e-07, maxiter=1500, trace=FALSE)

namc <- names(control)
if (!all(namc %in% names(control.default))) 
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
ctrl <- modifyList(control.default, control)

#
# method = reduced-rank ("rre") or minimal-polynomial ("mpe") extrapolation
# K = order of extrapolation scheme; K=2,3,4 are typical choices
# square = a logical variable indicating whether or not "squaring" is used

    method <- ctrl$method
    K <- ctrl$K
    square <- ctrl$square
    tol <- ctrl$tol
    maxiter <- ctrl$maxiter
    kr <- ctrl$kr
    trace <- ctrl$trace
 if (trace) cat("Cyclem-2 \n")

if ( K > length(par)) {
cat("K is too large.  Decrease it.  \n")
return()
}

iter <- 1

ptmp <- fixptfn(par,...)
res <- sqrt(c(crossprod(ptmp - par)))
par <- ptmp
conv <- TRUE
feval <- 1

while (feval < maxiter & res > tol) {
	extrap <- TRUE
	p <- matrix(NA, nrow=length(par), ncol=K+2)
	U <- matrix(NA, nrow=length(par), ncol=K+1)
	p[,1] <- par
	for (i in 2:(K+2)){
	p[,i] <- try(fixptfn(p[,i-1], ...), silent=TRUE)
	if (class(p[,i]) == "try-error" | any(is.nan(unlist(p[,i]))) | any(!is.numeric(p[,i]))) stop("Error in function evaluation")
	U[,i-1] <- p[, i] - p[, i-1]
	}
	feval <- feval + K + 1 
	sqKp1 <- sqrt(c(crossprod(U[, K+1])))

	if (method == "rre"){
	coef <- try(solve(qr(crossprod(U), LAPACK=TRUE, tol=1.e-14), rep(1,K+1)), silent=TRUE)
	if (class(coef) == "try-error" | any(is.nan(coef))) {
		extrap <- FALSE
	} else {
		if (abs(sum(coef)) < 1.e-07) extrap <- FALSE
		coef <- coef/sum(coef)
		}
	}
	if (method == "mpe"){
	coef <- try(solve(qr(U[,-(K+1)], LAPACK=TRUE, tol=1.e-14), -U[,K+1]), silent=TRUE)
	if (class(coef) == "try-error" | any(is.nan(coef))) {
		extrap <- FALSE
	} else {
		coef <- c(coef, 1)
		if (abs(sum(coef)) < 1.e-07) extrap <- FALSE
		coef <- coef/sum(coef)
		}
	}

	if (!extrap) {
		par <- p[, ncol(p)]
		res <- sqKp1
		iter <- iter + 1
 		if (trace) cat(" residual: ", res, "extrapolation: ", extrap, "\n")
		next 
	}

	if (!square) pnew <- c(p[,-(K+2)] %*% coef)
	if (square) {
		pnew <- rep(0, length(par))
		if (K > 1) {
		for (i in 1:(K-1)) {
			p <- try(cbind(p, fixptfn(p[,i+K+1],...)), silent=TRUE)
			if (class(p) == "try-error" | any(is.nan(unlist(p))) | any(!is.numeric(p))) stop("Error in function evaluation")
			}
		feval <- feval + K - 1
		}
		for (i in 0:K){
		for (j in 0:K){
			pnew <- pnew + coef[i+1]* coef[j+1] * p[,i+j+1]
		} }
		sqKp1 <- sqrt(c(crossprod(p[, 2*K] - p[, 2*K-1])))
	}

#  Perform one stabilization step between cycles
		ptemp <- try(fixptfn(pnew, ...), silent=TRUE)
		
      feval <- feval + 1
	if(class(ptemp)=="try-error" | any(is.nan(ptemp)) ){
		par <- p[, ncol(p)]
		res <- sqKp1
		extrap <- FALSE
	} else {
		res <- sqrt(c(crossprod(ptemp - pnew)))
		parnorm <- as.numeric(sqrt(crossprod(p[, ncol(p)])/length(par)))
		if (res <= (kr * (1 + parnorm) + sqKp1)) {
			par <- ptemp 
			extrap <- TRUE
			} else {
				par <- p[, ncol(p)]
				res <- sqKp1
				extrap <- FALSE
			}
	} 	

	iter <- iter+1
 	if (trace) cat(" residual: ", res, "extrapolation: ", extrap, "\n")

}  # Main loop complete
	if (feval >= maxiter) conv <- FALSE

return(list(par=par, value.objfn=NA, iter = iter, fpevals=feval, objfevals = 0, convergence=conv))

}
###################################################################

# The EM algorithm
#
fpiter <- function(par, fixptfn, objfn=NULL, control=list( ), ...){

control.default <- list(tol=1.e-07, maxiter=5000, trace=FALSE)
namc <- names(control)
if (!all(namc %in% names(control.default))) 
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
ctrl <- modifyList(control.default, control)

#
# method = reduced-rank ("rre") or minimal-polynomial ("mpe") extrapolation
# K = order of extrapolation scheme; K=2,3,4 are typical choices
# square = a logical variable indicating whether or not "squaring" is used

    tol <- ctrl$tol
    maxiter <- ctrl$maxiter
    trace <- ctrl$trace

if (trace) cat("fpiter \n")

iter <- 1
resid <- rep(NA,1)
objeval <- 0
conv <- FALSE

while (iter < maxiter) {

p.new <- fixptfn(par, ...)
res <- sqrt(crossprod(p.new - par))

if ( res < tol) {conv <- TRUE; break}

if (trace) {
     if (!is.null(objfn)) {cat("Iter: ", iter, "Objective fn: ",objfn(par, ...), "\n"); objeval <- objeval + 1}
     else cat("Iter: ", iter, "Residual: ",res, "\n")
	}
  par <- p.new
iter <- iter+1
}

loglik.best <-  if (!is.null(objfn)) objfn(par, ...) else NA

return(list(par=par, value.objfn=loglik.best, fpevals=iter, objfevals = objeval, convergence=conv))
}

###################################################################


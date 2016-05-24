constrOptim.nl <-
function (par, fn, gr=NULL, hin=NULL, hin.jac=NULL, heq=NULL, heq.jac=NULL, 
control.outer = list(), control.optim = list(), ...)  {

   if (is.null(heq) & is.null(hin)) stop("This is an unconstrained optimization problem - you should use `optim' \n")

control.outer.default <- list(mu0 = 0.01, sig0 = 10, eps = 1e-07,
       itmax = 50, method = "BFGS", trace = TRUE, NMinit = FALSE)

control.optim.default <- list(trace = 0, fnscale = 1, parscale = rep.int(1,
       length(par)), ndeps = rep.int(0.001, length(par)), maxit = 100L,
       abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1,
       beta = 0.5, gamma = 2, REPORT = 10, type = 1, lmm = 5,
       factr = 1e+07, pgtol = 0, tmax = 10, temp = 10)

outer.ctrl <- modifyList(control.outer.default, control.outer) 
optim.ctrl <- modifyList(control.optim.default, control.optim)

#	require(numDeriv, quietly=TRUE)
     if (is.null(gr)) gr <- function(par, ...) grad(func=fn, x=par, method= "simple", ...)

   if (is.null(hin)) {
    if (is.null(heq.jac) ) heq.jac <- function(par, ...) jacobian(func=heq, x=par, method= "simple", ...)
    ans <- augpen(par, fn, gr, heq=heq, heq.jac=heq.jac, control.outer=outer.ctrl, control.optim=optim.ctrl, ...)
  }  else if (is.null(heq)) {

    if (is.null(hin.jac)) hin.jac <- function(par, ...) jacobian(func=hin, x=par, method= "simple", ...)
    ans <- adpbar(par, fn, gr, hin=hin, hin.jac=hin.jac, control.outer=outer.ctrl, control.optim=optim.ctrl, ...)

  }   else  {
    if (is.null(heq.jac) ) heq.jac <- function(par, ...) jacobian(func=heq, x=par, method= "simple", ...)
    if (is.null(hin.jac)) hin.jac <- function(par, ...) jacobian(func=hin, x=par, method= "simple", ...)
	ans <- alabama(par, fn, gr, hin=hin, hin.jac=hin.jac, heq=heq, heq.jac=heq.jac, control.outer=outer.ctrl, control.optim=optim.ctrl, ...)
	}

    return(ans)
}

###################################################################
adpbar <-
function (theta, fn, gr=gr, hin=hin, hin.jac=hin.jac, control.outer = control.outer, control.optim=control.optim, ...)  {

mu <- control.outer$mu0
trace <- control.outer$trace
eps <- control.outer$eps
itmax <- control.outer$itmax
method <- control.outer$method
NMinit <- control.outer$NMinit

     if (!is.null(control.optim$fnscale) && control.optim$fnscale < 0) 
         mu <- -mu
 
    R <- function(theta, theta.old, ...) {
        gi <- hin(theta, ...)
        if (any(gi < 0)) return(NaN)
        gi.old <- hin(theta.old, ...)
	bar <- sum(gi.old * log(gi) - hin.jac(theta.old, ...) %*% theta)

        if (!is.finite(bar)) 
            bar <- -Inf
      fn(theta, ...) - mu * bar
    }

    dR <- function(theta, theta.old, ...) {
        gi <- hin(theta, ...)
	gi.old <- hin(theta.old, ...)
	hi <- hin.jac(theta.old, ...)         
        dbar <- colSums(hi* gi.old/gi - hi)
        gr(theta, ...) - mu * dbar
    }
   
    if (any(hin(theta, ...) <= 0)) 
        stop("initial value not feasible")

    obj <- fn(theta, ...)
    r <- R(theta, theta, ...)
	feval <- 0
	geval <- 0 

	  h0 <- hin(theta, ...)
  	mu <- mu * min(h0)
	if (abs(mu) < 1.e-10) mu <- 1.e-04 * sign(mu)

    for (i in 1:itmax) {
	if (trace) {
		cat("Min(hin): ", min(h0), "\n")
		cat("par: ", theta, "\n")
		cat("fval: ", obj, "\n")
	}
        obj.old <- obj
        r.old <- r
        theta.old <- theta

        fun <- function(theta, ...) {
            R(theta, theta.old, ...)
        }
        grad <- function(theta, ...) {
            dR(theta, theta.old, ...)
        }
	
       if ( NMinit & i == 1)  a <- optim(par=theta.old, fn=fun, gr=grad, control = control.optim, method = "Nelder-Mead", ...)
      else a <- optim(par=theta.old, fn=fun, gr=grad, control = control.optim, method = method, ...)
        r <- a$value

#  Here is "absolute" convergence criterion:
        if (is.finite(r) && is.finite(r.old) && abs(r - r.old) < eps) 
            break
        theta <- a$par
	  h0 <- hin(theta, ...)
	feval <- feval + a$counts[1]
	if (!NMinit | i > 1) geval <- geval + a$counts[2]
        obj <- fn(theta, ...)
        if (obj * sign(mu) > obj.old * sign(mu)) 
            break
    }

    if (i == itmax) {
        a$convergence <- 7
        a$message <- "Barrier algorithm ran out of iterations and did not converge"
    }
    if (mu > 0 && obj > obj.old) {
        a$convergence <- 11
        a$message <- paste("Objective function increased at outer iteration", 
            i)
    }
    if (mu < 0 && obj < obj.old) {
        a$convergence <- 11
        a$message <- paste("Objective function decreased at outer iteration", 
            i)
    }
    a$outer.iterations <- i
    a$barrier.value <- a$value
    a$value <- fn(a$par, ...)
    a$barrier.value <- a$barrier.value - a$value
    a$counts <- c(feval, geval)  # total number of fevals and gevals in inner iterations in BFGS
    a
}
###################################################################
augpen <-
function (theta, fn, gr=gr, heq=heq, heq.jac=heq.jac, control.outer = control.outer, control.optim=control.optim, ...)  {

mu <- control.outer$mu0
sig <- control.outer$sig0
trace <- control.outer$trace
eps <- control.outer$eps
itmax <- control.outer$itmax
method <- control.outer$method
NMinit <- control.outer$NMinit

   if (!is.null(control.optim$fnscale) && control.optim$fnscale < 0) {
		pfact <- -1
	} else pfact <- 1


	ans <- vector("list")
	feval <- 0
	geval <- 0 
  	k <- 0
	lam0 <- 0
	lam <- lam0
	i0 <- heq(theta, ...)
	Kprev <- max(abs(i0))
	sig0 <-  sig / Kprev 
	if (is.infinite(sig0)) sig0 <- 1
	sig <- sig0
	K <- Inf
	r <- fn(theta, ...)

	while (K > eps & k <= itmax) {   # Augmented Lagrangian loop for nonlinear "equality" constraintts 

	if (trace) {
		cat("Max(abs(heq)): ", Kprev, "\n")
		cat("par: ", theta, "\n")
		cat("fval: ", r, "\n")
	}

        fun <- function(theta, ...) {
		it <- heq(theta, ...)
            fn(theta, ...) - pfact * sum (lam * it) + pfact * sig/2 * sum(it * it)
        }
        grad <- function(theta, ...) {
		it <- heq(theta, ...)
		ij <- heq.jac(theta, ...)
            gr(theta, ...) - pfact * colSums(lam * ij) + pfact * sig * drop(t(ij) %*% it)
        }  
	
        if ( NMinit & k == 0 ) a <- optim(par=theta, fn=fun, gr=grad, control = control.optim, method = "Nelder-Mead", ...)
        else a <- optim(par=theta, fn=fun, gr=grad, control = control.optim, method = method, ...)

	  theta <- a$par
        r <- a$value
	  i0 <- heq(theta, ...)
	  K <- max(abs(i0)) 
        feval <- feval + a$counts[1]
	  if (!NMinit | k > 0) geval <- geval + a$counts[2]
	  k <- k + 1

		if( K <= Kprev/4) {
			lam <- lam - i0 * sig
			Kprev <- K
		} else sig <- 10 * sig
	}  # Augmented Lagrangian loop completed

     if (k >= itmax) {
        a$convergence <- 7
        a$message <- "Augmented Lagrangian algorithm ran out of iterations and did not converge"
    }
  
    ans$par <- theta
    ans$value <- fn(a$par, ...)
    ans$iterations <- k
    ans$lambda <- lam
    ans$penalty <- r - ans$value
    ans$counts <- c(feval, geval)  # total number of fevals and gevals in inner iterations in BFGS

    ans
}

###################################################################
alabama <-
function (theta, fn, gr=gr, hin=hin, hin.jac=hin.jac, heq=heq, heq.jac=heq.jac, control.outer = control.outer, control.optim=control.optim, ...) {

mu <- control.outer$mu0
sig <- control.outer$sig0
trace <- control.outer$trace
eps <- control.outer$eps
itmax <- control.outer$itmax
method <- control.outer$method
NMinit <- control.outer$NMinit

	if (is.null(heq) & is.null(hin)) stop("This is an unconstrained optimization problem - use `optim' \n")

    if (!is.null(control.optim$fnscale) && control.optim$fnscale < 0) {
		mu <- -mu
		pfact <- -1
	} else pfact <- 1

	alpha <- 0.5
	beta <- 1

    R <- function(theta, theta.old, ...) {
        gi <- hin(theta, ...)
        if (any(gi < 0)) return(NaN)
        gi.old <- hin(theta.old, ...)
	  hjac <- hin.jac(theta.old, ...)
	bar <- sum(gi.old * log(gi) - hjac %*% theta)

        if (!is.finite(bar)) 
            bar <- -Inf
      fn(theta, ...) - mu * bar
    }

    dR <- function(theta, theta.old, ...) {
        gi <- hin(theta, ...)
	gi.old <- hin(theta.old, ...)
	hjac <- hin.jac(theta.old, ...)         
        dbar <- colSums(hjac* gi.old/gi - hjac)
        gr(theta, ...) - mu * dbar
    }
   
	h0 <- hin(theta, ...)
     if (any(h0 <= 0)) 
        stop("initial value violates inequality constraints")

    obj <- fn(theta, ...)
    r <- R(theta, theta, ...)
	feval <- 0
	geval <- 0 
	i0 <- heq(theta, ...)
	Kprev <- max(abs(i0))
	lam0 <- 0
	sig0 <-  sig / Kprev 
	if (is.infinite(sig0)) sig0 <- 1
	lam <- lam0
	sig <- sig0
	mu <- mu * min(h0)
	if (abs(mu) < 1.e-10) mu <- 1.e-04 * sign(mu)

	K <- Inf
	if (trace) cat("Min(hin): ", min(h0), "Max(abs(heq)): ", Kprev, "\n")

    for (i in 1:itmax) {  # Adaptive Barrier MM loop for nonlinear "inequality" constraintts 
	
	if (trace) {
		cat("Outer iteration: ", i, "\n")
		cat("Min(hin): ", min(h0), "Max(abs(heq)): ", Kprev, "\n")
		cat("par: ", signif(theta,6), "\n")
		cat("fval =  ", signif(obj,4), "\n \n")
	}

        obj.old <- obj
        r.old <- r
        theta.old <- theta

        fun <- function(theta, ...) {
		it <- heq(theta, ...)
            R(theta, theta.old, ...) - pfact * sum (lam * it) + pfact * sig/2 * sum(it * it)
        }

        grad <- function(theta, ...) {
		it <- heq(theta, ...)
		ij <- heq.jac(theta, ...)
            dR(theta, theta.old, ...) - pfact * colSums(lam * ij) + pfact * sig * drop(t(ij) %*% it)
        }  

	if(sig > 1e05) control.optim$reltol <- 1.e-10

        if ( NMinit & i == 1)  a <- optim(par=theta.old, fn=fun, gr=grad, control = control.optim, method = "Nelder-Mead", ...)
       else a <- optim(par=theta.old, fn=fun, gr=grad, control = control.optim, method = method, ...)

	theta <- a$par
        r <- a$value

	  h0 <- hin(theta, ...)
	  i0 <- heq(theta, ...)
	  K <- max(abs(i0))
         feval <- feval + a$counts[1]
	if (!NMinit | i > 1) geval <- geval + a$counts[2]

		if( K <= Kprev/4) {
			lam <- lam - i0 * sig
			Kprev <- K
		} else sig <- 10 * sig

        obj <- fn(theta, ...)

 #  Here is "absolute" convergence criterion:
	pconv <- max(abs(theta - theta.old))
        if ((is.finite(r) && is.finite(r.old) && abs(r - r.old) < eps && K < eps) | pconv < 1.e-12) {
#        if (is.finite(obj) && is.finite(obj.old) && abs(obj - obj.old) < eps && K < eps) {
		theta.old <- theta
		atemp <- optim(par=theta, fn=fun, gr=grad, control = control.optim, method = "BFGS", hessian=TRUE, ...)
		a$hessian <- atemp$hess 
		break
		}
		
}

    if (i == itmax) {
        a$convergence <- 7
        a$message <- "ALABaMA ran out of iterations and did not converge"
	}

    a$outer.iterations <- i
    a$lambda <- lam
    a$sigma <- sig
    a$barrier.value <- a$value
    a$value <- fn(a$par, ...)
    a$barrier.value <- a$barrier.value - a$value
    a$K <- K
    a$counts <- c(feval, geval)  # total number of fevals and gevals in inner iterations in BFGS
    a
}





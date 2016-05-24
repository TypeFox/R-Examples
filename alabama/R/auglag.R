
#####################################################################################
auglag <- function (par, fn, gr, hin, hin.jac, heq, heq.jac, 
	control.outer = list(), control.optim = list(), 
    ...) 
{
   if (missing(heq) & missing(hin)) stop("This is an unconstrained optimization problem - you should use `optim' \n")

control.outer.default <- list(lam0 = 10, sig0 = 100, eps = 1e-07,
       itmax = 50, method = "BFGS", trace = TRUE, NMinit = FALSE, ilack.max=6,
i.scale = 1, e.scale = 1, kkt2.check=TRUE)

control.optim.default <- list(trace = 0, fnscale = 1, parscale = rep.int(1,
       length(par)), ndeps = rep.int(0.001, length(par)), maxit = 100L,
       abstol = -Inf, reltol = sqrt(.Machine$double.eps), alpha = 1,
       beta = 0.5, gamma = 2, REPORT = 10, type = 1, lmm = 5,
       factr = 1e+07, pgtol = 0, tmax = 10, temp = 10)

control.nlminb.default <- list(eval.max=200,iter.max=150,trace=0,abs.tol=0,rel.tol=1.e-10,x.tol=1.5e-08,xf.tol=2.2e-14,step.min=1,step.max=1,
sing.tol=1.e-10)

control.outer <- modifyList(control.outer.default, control.outer) 

control.inner <- if (control.outer$method != "nlminb") modifyList(control.optim.default,control.optim) else if (control.outer$method=="nlminb") modifyList(control.nlminb.default, control.optim)

if (control.outer$method == "nlminb" & control.outer$NMinit) {
warning("NM initialization can't be used with `nlminb'; it will be ignored")
control.outer$NMinit <- FALSE
}

e.scale <- control.outer$e.scale
i.scale <- control.outer$i.scale

#	require(numDeriv, quietly=TRUE)
     if (missing(gr)) gr <- function(par, ...) grad(func=fn, x=par, method= "simple", ...) 

   if (missing(hin)) {
    heq.scaled <- function(par, ...) { heq(par, ...) / e.scale}
    heq.jac.scaled <- if (missing(heq.jac)) function(par, ...) jacobian(func=heq.scaled, x=par, method= "simple", ...) 
	else function(par, ...) heq.jac(par, ...) / e.scale

	ans <- auglag1(par, fn, gr, heq.scaled, heq.jac.scaled,  
	control.outer = control.outer, control.optim = control.inner, ...)  
	}  else if (missing(heq)) {
	hin.scaled <- function(par, ...) { hin(par, ...) / i.scale}
      hin.jac.scaled <-  if (missing(hin.jac)) function(par, ...) jacobian(func=hin.scaled, x=par, method= "simple", ...) 
	else function(par, ...) hin.jac(par, ...) / i.scale

	ans <- auglag2(par, fn, gr, hin.scaled, hin.jac.scaled,  
	control.outer = control.outer, control.optim = control.inner, ...) 
  }   else  {
    	heq.scaled <- function(par, ...) { heq(par, ...) / e.scale}
	hin.scaled <- function(par, ...) { hin(par, ...) / i.scale}
     	heq.jac.scaled <- if (missing(heq.jac) )function(par, ...) 			jacobian(func=heq.scaled, x=par, method= "simple", ...) 
	else function(par, ...) heq.jac(par, ...) / e.scale
    hin.jac.scaled <- if (missing(hin.jac)) function(par, ...) jacobian(func=hin.scaled, x=par, method= "simple", ...) 
	else function(par, ...) hin.jac(par, ...) / i.scale

	ans <- auglag3(par, fn=fn, gr=gr, hin=hin.scaled, hin.jac=hin.jac.scaled, heq=heq.scaled, heq.jac=heq.jac.scaled,  
	control.optim = control.inner, control.outer = control.outer, ...) 
	}

if (!missing(hin)) ans$ineq <- ans$ineq * i.scale
if (!missing(heq)) ans$equal <- ans$equal * e.scale

return(ans)
}

##################################################################
auglag1 <- function (par, fn, gr = NULL, 
	heq = NULL, heq.jac = NULL, 
	control.outer = list(), control.optim = list(), 
    ...) 
{

    sig <- control.outer$sig0
    lam0 <- control.outer$lam0
    trace <- control.outer$trace
    eps <- control.outer$eps
    itmax <- control.outer$itmax
    ilack.max <- control.outer$ilack.max
    method <- control.outer$method
    NMinit <- control.outer$NMinit
    kkt2.check <- control.outer$kkt2.check
    pfact <- if (!is.null(control.optim$fnscale) && control.optim$fnscale < 
        0) -1 else 1

        fun <- function(par, ...) {
		d0 <- heq(par, ...)
            fn(par, ...) - pfact * sum(lam * d0) + 
                pfact * sig/2 * sum(d0 * d0)
        }

        gradient <- function(par, ...) {
		d0 <- heq(par, ...)
            ij <- heq.jac(par, ...)
            gr(par, ...) - pfact * colSums(lam * ij) + 
	pfact * sig * drop(t(ij) %*% d0)
        }

    d0 <- heq(par, ...)
    lam <- rep(lam0, length(d0))
    		dmax <- max(abs(d0))

    obj <- fn(par, ...)
	r <- obj
    feval <- 0
    geval <- 0
    ilack <- 0
    Kprev <- dmax
    sig0 <- sig/Kprev
    if (is.infinite(sig0)) 
        sig0 <- 1
    sig <- sig0

    K <- Inf
    if (trace) cat("Max(abs(heq)): ", max(abs(d0)), "\n")
    for (i in 1:itmax) {
        if (trace) {
            cat("Outer iteration: ", i, "\n")
            cat("Max(abs(heq)): ", max(abs(d0)), "\n")
            cat("par: ", signif(par, 6), "\n")
            cat("fval =  ", signif(obj, 4), "\n \n")
        }
        par.old <- par
        obj.old <- obj
	  r.old <- r
        if (NMinit & i == 1) 
            a <- optim(par = par, fn = fun, 
                control = control.optim, method = "Nelder-Mead", ...)
        else if (method != "nlminb") {
        if (sig > 1e+05) 
            control.optim$reltol <- 1e-10
        a <- optim(par = par, fn = fun, gr=gradient, control = control.optim, method = method, ...)
        r <- a$value
        feval <- feval + a$counts[1]
        geval <- geval + a$counts[2]
        } else if (method == "nlminb") {
        a <- nlminb(start = par, objective = fun, gradient=gradient, control=control.optim, ...)
        r <- a$objective
        feval <- feval + a$evaluations[1]
        geval <- geval + a$evaluations[2]
        }
        par <- a$par
        d0 <- heq(par, ...)
        K <- max(abs(d0))
        if (K <= Kprev/4) {
           lam <- lam - d0 * sig
            Kprev <- K
        }
        else sig <- 10 * sig
        obj <- fn(par, ...)

        pconv <- max(abs(par - par.old))
	if (pconv < eps) {
		ilack <- ilack + 1
	} else ilack <- 0 

        if ((is.finite(r) && is.finite(r.old) && abs(r - r.old) < 
            eps && K < eps) | ilack >= ilack.max) break
}

    if (i == itmax) {
        a$convergence <- 7
        a$message <- "ALABaMA ran out of iterations and did not converge"
    } 
	else if (K > eps) {
            a$convergence <- 9
            a$message <- "Convergence due to lack of progress in parameter updates"
    }
    if (method=="nlminb") {
       a$evaluations <- NULL
       a$objective <- NULL
       a$iterations <- NULL
    }
    a$outer.iterations <- i
    a$lambda <- lam
    a$sigma <- sig
    a$value <- fn(a$par, ...)
    a$gradient <- gradient(a$par, ...)
    a$ineq <- NA
    a$equal <- heq(a$par, ...)
    a$counts <- c(feval, geval)
    a$kkt1 <- max(abs(a$gradient)) <= 0.01 * (1 + abs(a$value)) 
     a$kkt2 <- NA
    if (kkt2.check) {
    a$hessian <- jacobian(x=a$par, func=gradient, ...)
    evs <- eigen(a$hessian, symmetric=TRUE, only.values=TRUE)$value    
#    a$kkt2 <- any(evs * control.optim$fnscale > 0) 
    a$kkt2 <- if (any(abs(Im(evs)) > 1.e-14)) FALSE else if (all(Re(evs) * control.optim$fnscale > 0)) TRUE else FALSE  
    }
    a
}

##################################################################
auglag2 <- function (par, fn, gr = NULL, 
	hin = NULL, hin.jac = NULL, 
	control.outer = list(), control.optim = list(), 
    ...) 
{
    sig <- control.outer$sig0
    lam0 <- control.outer$lam0
    trace <- control.outer$trace
    eps <- control.outer$eps
    itmax <- control.outer$itmax
    ilack.max <- control.outer$ilack.max
    method <- control.outer$method
    NMinit <- control.outer$NMinit
    kkt2.check <- control.outer$kkt2.check
    pfact <- if (!is.null(control.optim$fnscale) && control.optim$fnscale < 
        0) -1 else 1

        fun <- function(par, ...) {
		h0 <- hin(par, ...)
    		d0 <- h0
    		inactive <- (1:length(h0))[(h0 > lam[1:length(h0)]/sig)]
		d0[inactive] <- lam[inactive] / sig
            fn(par, ...) - pfact * sum(lam * d0) + 
                pfact * sig/2 * sum(d0 * d0)
        }

        gradient <- function(par, ...) {
		h0 <- hin(par, ...)
		d0 <- h0
	    	active <- (1:length(h0))[(h0 <= lam[1:length(h0)]/sig)]
            ij <- hin.jac(par, ...)[active, , drop=FALSE]
		gr(par, ...) - pfact * colSums(lam[active] * 
                ij) + pfact * sig * drop(crossprod(ij, d0[active]))
        }

    	h0 <- hin(par, ...)
	d0 <- h0
    	lam <- rep(lam0, length(d0))
	    	inactive <- (1:length(h0))[(h0 > lam[1:length(h0)]/sig)]
		d0[inactive] <- lam[inactive] / sig
    		dmax <- max(abs(d0))

    obj <- fn(par, ...)
	r <- obj
    feval <- 0
    geval <- 0
    ilack <- 0
    Kprev <- dmax
    sig0 <- sig/Kprev
    if (is.infinite(sig0)) 
        sig0 <- 1
    sig <- sig0

    K <- Inf
    if (trace) cat("Min(hin): ", min(h0), "\n")
    for (i in 1:itmax) {
        if (trace) {
            cat("Outer iteration: ", i, "\n")
            cat("Min(hin): ", min(h0), "\n")
            cat("par: ", signif(par, 6), "\n")
            cat("fval =  ", signif(obj, 4), "\n \n")
        }
        par.old <- par
        obj.old <- obj
	  r.old <- r
        if (NMinit & i == 1) 
            a <- optim(par = par, fn = fun, 
                control = control.optim, method = "Nelder-Mead", ...)
        else if (method != "nlminb") {
        if (sig > 1e+05) 
            control.optim$reltol <- 1e-10
        a <- optim(par = par, fn = fun, gr=gradient, control = control.optim, method = method, ...)
        r <- a$value
        feval <- feval + a$counts[1]
        geval <- geval + a$counts[2]
        } else if (method == "nlminb") {
        a <- nlminb(start = par, objective = fun, gradient=gradient, control=control.optim, ...)
        r <- a$objective
        feval <- feval + a$evaluations[1]
        geval <- geval + a$evaluations[2]
        }
        par <- a$par
    	h0 <- hin(par, ...)
	d0 <- h0
	    	inactive <- (1:length(h0))[(h0 > lam[1:length(h0)]/sig)]
		d0[inactive] <- lam[inactive] / sig
        K <- max(abs(d0))
        if (K <= Kprev/4) {
           lam <- lam - d0 * sig
            Kprev <- K
        }
        else sig <- 10 * sig
        obj <- fn(par, ...)

        pconv <- max(abs(par - par.old))
	if (pconv < eps) {
		ilack <- ilack + 1
	} else ilack <- 0 

        if ((is.finite(r) && is.finite(r.old) && abs(r - r.old) < 
            eps && K < eps) | ilack >= ilack.max) break
}

    if (i == itmax) {
        a$convergence <- 7
        a$message <- "ALABaMA ran out of iterations and did not converge"
    } 
	else if (K > eps) {
            a$convergence <- 9
            a$message <- "Convergence due to lack of progress in parameter updates"
    }
    if (method=="nlminb") {
       a$evaluations <- NULL
       a$objective <- NULL
       a$iterations <- NULL
    }
    a$outer.iterations <- i
    a$lambda <- lam
    a$sigma <- sig
    a$value <- fn(a$par, ...)
    a$gradient <- gradient(a$par, ...)
    a$ineq <- hin(a$par, ...) 
    a$equal <- NA
    a$counts <- c(feval, geval)
    a$kkt1 <- max(abs(a$gradient)) <= 0.01 * (1 + abs(a$value)) 
     a$kkt2 <- NA
    if (kkt2.check) {
    a$hessian <- jacobian(x=a$par, func=gradient, ...)
    evs <- eigen(a$hessian, symmetric=TRUE, only.values=TRUE)$value    
#    a$kkt2 <- any(evs * control.optim$fnscale > 0) 
    a$kkt2 <- if (any(abs(Im(evs)) > 1.e-14)) FALSE else if (all(Re(evs) * control.optim$fnscale > 0)) TRUE else FALSE  
    }
    a
}

##################################################################
auglag3 <- function (par, fn, gr = NULL, 
	hin = NULL, hin.jac = NULL, heq = NULL, heq.jac = NULL, 
	control.outer = list(), control.optim = list(), 
    ...) 
{

    sig <- control.outer$sig0
    lam0 <- control.outer$lam0
    trace <- control.outer$trace
    eps <- control.outer$eps
    itmax <- control.outer$itmax
    ilack.max <- control.outer$ilack.max
    method <- control.outer$method
    NMinit <- control.outer$NMinit
    kkt2.check <- control.outer$kkt2.check
    pfact <- if (!is.null(control.optim$fnscale) && control.optim$fnscale < 
        0) -1 else 1

        fun <- function(par, ...) {
	h0 <- hin(par, ...)
	i0 <- heq(par, ...)
    	d0 <- c(h0, i0)
	    	active <- (1:length(h0))[(h0 <= lam[1:length(h0)]/sig)]
	    	inactive <- (1:length(h0))[(h0 > lam[1:length(h0)]/sig)]
		d0[active] <- h0[active]
		d0[inactive] <- lam[inactive] / sig
            fn(par, ...) - pfact * sum(lam * d0) + 
                pfact * sig/2 * sum(d0 * d0)
        }

        gradient <- function(par, ...) {
	h0 <- hin(par, ...)
	i0 <- heq(par, ...)
    	d0 <- c(h0, i0)
	    	active <- (1:length(h0))[(h0 <= lam[1:length(h0)]/sig)]
	    	inactive <- (1:length(h0))[(h0 > lam[1:length(h0)]/sig)]
            ij <- rbind(hin.jac(par, ...)[active, , drop=FALSE], heq.jac(par, ...))
	if (length(inactive) > 0) 
		gr(par, ...) - pfact * colSums(lam[-inactive] * ij) + 
pfact * sig * drop(t(ij) %*% d0[-inactive])
	else 
			gr(par, ...) - pfact * colSums(lam * ij) + 
			pfact * sig * drop(t(ij) %*% d0)
        }

    h0 <- hin(par, ...)
    i0 <- heq(par, ...)
    d0 <- c(h0, i0)

    lam <- rep(lam0, length(d0))
	    	active <- (1:length(h0))[(h0 <= lam[1:length(h0)]/sig)]
	    	inactive <- (1:length(h0))[(h0 > lam[1:length(h0)]/sig)]
		d0[active] <- h0[active]
		d0[inactive] <- lam[inactive] / sig
    		dmax <- max(abs(d0))

    obj <- fn(par, ...)
	r <- obj
    feval <- 0
    geval <- 0
    ilack <- 0
    Kprev <- dmax
    sig0 <- sig/Kprev
    if (is.infinite(sig0)) 
        sig0 <- 1
    sig <- sig0

    K <- Inf
    if (trace) cat("Min(hin): ", min(h0), "Max(abs(heq)): ", max(abs(i0)), "\n")
    for (i in 1:itmax) {
        if (trace) {
            cat("Outer iteration: ", i, "\n")
            cat("Min(hin): ", min(h0), "Max(abs(heq)): ", max(abs(i0)), 
                "\n")
            cat("par: ", signif(par, 6), "\n")
            cat("fval =  ", signif(obj, 4), "\n \n")
        }
        par.old <- par
        obj.old <- obj
	  r.old <- r
        if (NMinit && (i == 1)) 
            a <- optim(par = par, fn = fun, 
                control = control.optim, method = "Nelder-Mead", ...)
        else if (method != "nlminb") {
        if (sig > 1e+05) 
            control.optim$reltol <- 1e-10
        a <- optim(par = par, fn = fun, gr=gradient, control = control.optim, method = method, ...)
        r <- a$value
        feval <- feval + a$counts[1]
        geval <- geval + a$counts[2]
        } else if (method == "nlminb") {
        a <- nlminb(start = par, objective = fun, gradient=gradient, control=control.optim, ...)
        r <- a$objective
        feval <- feval + a$evaluations[1]
        geval <- geval + a$evaluations[2]
        }
        par <- a$par
        h0 <- hin(par, ...)
        i0 <- heq(par, ...)
		d0 <- c(h0, i0)
	    	active <- (1:length(h0))[(h0 <= lam[1:length(h0)]/sig)]
	    	inactive <- (1:length(h0))[(h0 > lam[1:length(h0)]/sig)]
		if (length(active) > 0) d0[active] <- h0[active]
		if (length(inactive) > 0) d0[inactive] <- lam[inactive]/sig
        K <- max(abs(d0))
        if (K <= Kprev/4) {
           lam <- lam - d0 * sig
            Kprev <- K
        }
        else sig <- 10 * sig
        obj <- fn(par, ...)

        pconv <- max(abs(par - par.old))
	if (pconv < eps) {
		ilack <- ilack + 1
	} else ilack <- 0 

        if ((is.finite(r) && is.finite(r.old) && abs(r - r.old) < 
            eps && K < eps) | ilack >= ilack.max) break
}

    if (i == itmax) {
        a$convergence <- 7
        a$message <- "ALABaMA ran out of iterations and did not converge"
    } 
	else if (K > eps) {
            a$convergence <- 9
            a$message <- "Convergence due to lack of progress in parameter updates"
    }
    if (method=="nlminb") {
       a$evaluations <- NULL
       a$objective <- NULL
       a$iterations <- NULL
    }
    a$outer.iterations <- i
    a$lambda <- lam
    a$sigma <- sig
    a$value <- fn(a$par, ...)
    a$gradient <- gradient(a$par, ...)
    a$ineq <- hin(a$par, ...) 
    a$equal <- heq(a$par, ...)
    a$counts <- c(feval, geval)
    a$kkt1 <- max(abs(a$gradient)) <= 0.01 * (1 + abs(a$value)) 
     a$kkt2 <- NA
    if (kkt2.check) {
    a$hessian <- jacobian(x=a$par, func=gradient, ...)
    evs <- eigen(a$hessian, symmetric=TRUE, only.values=TRUE)$value    
#    a$kkt2 <- any(evs * control.optim$fnscale > 0) 
    a$kkt2 <- if (any(abs(Im(evs)) > 1.e-14)) FALSE else if (all(Re(evs) * control.optim$fnscale > 0)) TRUE else FALSE  
    }
    a
}



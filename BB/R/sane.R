sane <- function (par, fn, method = 2, control = list(), 
                  quiet=FALSE, alertConvergence=TRUE, ...)  {
    ctrl <- list(maxit = 1500, M = 10, tol = 1e-07, trace = TRUE, 
        triter = 10, quiet=FALSE, noimp = 100, NM=FALSE, BFGS=FALSE)
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    M <- ctrl$M
    maxit <- ctrl$maxit
    tol <- ctrl$tol
    trace <- ctrl$trace
    triter <- ctrl$triter
    noimp <- ctrl$noimp
    NM <- ctrl$NM
    BFGS <- ctrl$BFGS
    fargs <- list(...)

    lineSearch <- function(x, fn, F, fval, dg, M, lastfv, sgn, 
        lambda, fcnt, bl, fargs) {
        maxbl <- 100
        gamma <- 1e-04
        sigma1 <- 0.1
        sigma2 <- 0.5
        cbl <- 0
        fmax <- max(lastfv)
        gpd <- -2 * abs(dg)
        while (cbl < maxbl) {
            xnew <- x + lambda * sgn * F
            Fnew <- try(do.call(fn, append(list(xnew), fargs)))
            fcnt = fcnt + 1
            if (class(Fnew) == "try-error" || any(is.nan(Fnew))) 
                return(list(xnew = NA, Fnew = NA, fcnt = fcnt, 
                  bl = bl, lsflag = 1, fune = NA))
            else fune <- sum(Fnew * Fnew)
            if (fune <= (fmax + lambda * gpd * gamma)) {
                if (cbl >= 1) 
                  bl <- bl + 1
                return(list(xnew = xnew, Fnew = Fnew, fcnt = fcnt, 
                  lambda = lambda, bl = bl, lsflag = 0, fune = fune))
            }
            else {
                lamc <- -(gpd * lambda^2)/(2 * (fune - fval - 
                  lambda * gpd))
                c1 <- sigma1 * lambda
                c2 <- sigma2 * lambda
                if (lamc < c1) 
                  lambda <- c1
                else if (lamc > c2) 
                  lambda <- c2
                else lambda <- lamc
                cbl <- cbl + 1
            }
        }
        return(list(xnew = NA, Fnew = NA, fcnt = fcnt, lambda = NA, 
            bl = bl, lsflag = 2, fune = NA))
    }

    n <- length(par)
    fcnt <- iter <- bl <- 0
    alfa <- 1
    eps <- 1e-10
    h <- 1.e-07
    lastfv <- rep(0, M)

	U <- function(x, ...) drop(crossprod(fn(x, ...)))
##  We do initial Nelder-Mead start-up
	if (NM) {
		res <- try(optim(par=par, fn=U, method="Nelder-Mead", control=list(maxit=100), ...), silent=TRUE)
		if (class(res) == "try-error") { 
			cat(res)
			stop("\nFailure in Nelder-Mead Start.  Try another starting value \n")
			}
		else if (any(is.nan(res$par))) 
			stop("Failure in Nelder-Mead Start (NaN value).  Try another starting value \n")
		par <- res$par
		fcnt <- as.numeric(res$counts[1])
	}

    F <- try(fn(par, ...))
	fcnt <- fcnt + 1
    if (class(F) == "try-error") 
        stop(" Failure in initial functional evaluation.")
    else if (!is.numeric(F) || !is.vector(F)) 
        stop("Function must return a vector numeric value.")
    else if (any(is.nan(F), is.infinite(F), is.na(F))) 
        stop(" Failure in initial functional evaluation.")
    else if (length(F) == 1) if (!quiet)
        warning("Function returns a scalar. Function BBoptim or spg is better.")

    F0 <- normF <- sqrt(sum(F * F))
    if (trace) 
        cat("Iteration: ", 0, " ||F(x0)||: ", F0/sqrt(n), "\n")
    pbest <- par
    normF.best <- normF
    lastfv[1] <- normF^2
    flag <- 0
    knoimp <- 0

    while (normF/sqrt(n) > tol & iter <= maxit) {
        Fa <- try(fn(par + h * F, ...))
        fcnt <- fcnt + 1
        if (class(Fa) == "try-error" || any(is.nan(Fa))) {
            flag <- 1
            break
        }
        dg <- (sum(F * Fa) - normF^2)/h
        if (abs(dg/normF^2) < eps | is.nan(dg) | is.infinite(dg)) {
            flag <- 3
            break
        }
        if ((alfa <= eps) | (alfa >= 1/eps)) 
            alfa <- if (normF > 1) 
                1

            else if (normF >= 1e-05 & normF <= 1) 
                normF
            else if (normF < 1e-05) 
                1e-05
        sgn <- if (dg > 0) 
            -1
        else 1

######## change made on Aug 29, 2008
## Steplength for first iteration is scaled properly
##

	if (iter==0) {
	alfa <- max(normF, 1)
	alfa1 <- alfa2 <- alfa
	}
	
	temp <- alfa2
	alfa2 <- alfa
	if (normF <= 0.01) alfa <- alfa1
	alfa1 <- temp

        lambda <- 1/alfa
        ls.ret <- lineSearch(x = par, fn = fn, F = F, fval = normF^2, 
            dg = dg, M = M, lastfv = lastfv, sgn, lambda, fcnt, 
            bl, fargs)
        fcnt <- ls.ret$fcnt
        bl <- ls.ret$bl
        flag <- ls.ret$lsflag

       if (flag > 0) 
            break

        lambda <- ls.ret$lambda
        Fnew <- ls.ret$Fnew
        pnew <- ls.ret$xnew
        fune <- ls.ret$fune
        alfa <- if (method == 1) 
            sum(F * (F - Fnew))/(lambda * sum(F * F))
        else if (method == 2) 
            sum((F - Fnew)^2)/(lambda * sum(F * (F - Fnew)))
        else if (method == 3) 
            sqrt(sum((F - Fnew)^2)/(lambda^2 * sum(F * F)))

         if (is.nan(alfa))  
		alfa <- eps

       par <- pnew
        F <- Fnew
        fun <- fune
        normF <- sqrt(fun)

       if (normF < normF.best) {
            pbest <- par
            normF.best <- normF
		knoimp <- 0
        } else knoimp <- knoimp + 1

	if (knoimp == noimp) {
		flag <- 4
		break
		}

        iter <- iter + 1
        lastfv[1 + iter%%M] <- fun
        if (trace && (iter%%triter == 0)) 
            cat("\n iteration: ", iter, " ||F(xn)|| =  ", normF, 
                "\n")
    }
    if (flag == 0) {
        if (normF.best/sqrt(n) <= tol) 
            conv <- list(type = 0, message = "Successful convergence")
        if (iter > maxit) 
            conv <- list(type = 1, message = "Maximum number of iterations exceeded")
    }
    else if (flag == 1) 
            conv <- list(type = 2, message = "Error in function evaluation")
    else if (flag == 2) 
            conv <- list(type = 3, message = "Maximum limit on steplength reductions exceeded")
    else if (flag == 3) 
            conv <- list(type = 4, message = "Anomalous iteration")
    else if (flag == 4) 
            conv <- list(type = 5, message = "Lack of improvement in objective function")

##  We do final "optim" iterations using "L-BFGS-B" when type=4 or 5
	if (BFGS & (conv$type==4 | conv$type==5) ) {
	if (!quiet) cat("Calling `L-BFGS-B' in `optim' \n")
	res <- try(optim(par=pbest, fn=U, method="L-BFGS-B", 
	               control=list(pgtol=1.e-08, factr=1000, maxit=200), ...), 
	           silent=TRUE)
        if (!inherits(res, "try-error") && !any(is.nan(res$par)) ) {
      		normF.new <- sqrt(res$value)
 		if (normF.new < normF.best) {
			normF.best <- normF.new
			pbest <- res$par
			}
 		}
	fcnt <- fcnt + as.numeric(res$counts[1])

        if (normF.best/sqrt(length(par)) <= tol) 
        	conv <- list(type = 0, message = "Successful convergence")
	}
	
    if(alertConvergence && ( 0 != conv$type))
          warning("Unsuccessful convergence.")

    return(list(par = pbest, residual = normF.best/sqrt(length(par)), 
        fn.reduction = F0 - normF.best, feval = fcnt, iter = iter, 
        convergence = conv$type, message = conv$message))
}




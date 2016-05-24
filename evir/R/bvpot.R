"gpdbiv" <- 
function(data1 = NA, data2 = NA, u1 = NA, u2 = NA, ne1 = NA,
    ne2 = NA, global = FALSE, method = "BFGS", ...)
{
    data1 <- as.numeric(data1)
    data2 <- as.numeric(data2)
    
    Zfunc <- function(y, u, lambda, xi, sigma)
        (lambda^-1) * (1 + (xi * pmax((y - u), 0))/sigma)^(1/xi)
    Kfunc <- function(y, u, lambda, xi, sigma)
        -lambda^(-xi) * (sigma^-1) * (Zfunc(y, u, lambda, xi, sigma))^(1 - xi)
    Vfunc <- function(x, y, alpha)
        (x^(-1/alpha) + y^(-1/alpha))^alpha
    Vfunc1 <- function(x, y, alpha)
        -x^(-(1/alpha) - 1) * (x^(-1/alpha) + y^(-1/alpha))^(alpha - 1)
    Vfunc2 <- function(x, y, alpha)
        -(alpha - 1) * (alpha^-1) * (x * y)^(-(1/alpha) - 1) *
          (x^(-1/alpha) + y^(-1/alpha))^(alpha - 2)
    fun <- list(Z = Zfunc, K = Kfunc, V = Vfunc, V1 = Vfunc1, V2 = Vfunc2) 
   
    if(is.na(ne1) && is.na(u1))
        stop(paste("Enter either a threshold or",
                   "the number of upper extremes for margin 1"))
    if(!is.na(ne1) && !is.na(u1))
        stop(paste("Enter EITHER a threshold or",
                   "the number of upper extremes for margin 1"))
    if(is.na(ne2) && is.na(u2))
        stop(paste("Enter either a threshold or",
                   "the number of upper extremes for margin 2"))
    if(!is.na(ne2) && !is.na(u2))
        stop(paste("Enter EITHER a threshold or",
                   "the number of upper extremes for margin 2"))

    out1 <- gpd(data1, threshold = u1, nextremes = ne1)
    par.ests1 <- out1$par.ests
    par.ses1 <- out1$par.ses
    
    out2 <- gpd(data2, threshold = u2, nextremes = ne2)
    par.ests2 <- out2$par.ests
    par.ses2 <- out2$par.ses

    uu <- c(out1$threshold, out2$threshold)
    ne <- c(out1$n.exceed, out2$n.exceed)
    mpar <- c(par.ests1, par.ests2) 
    
    delta1 <- as.numeric(data1 > uu[1])
    delta2 <- as.numeric(data2 > uu[2])
    lambda1 <- sum(delta1)/length(data1)
    lambda2 <- sum(delta2)/length(data2)

    theta <- 0.8
    if(global) {
        theta <- c(theta, mpar)
        mpar <- NULL
    }
	
    negloglik <- function(theta, data1, data2, uu, delta1, delta2,
        lambda1, lambda2, mpar, fun)
    {
      	alpha <- theta[1]
	if(is.null(mpar)) {
            xi1 <- theta[2] ; sigma1 <- theta[3]
	    xi2 <- theta[4] ; sigma2 <- theta[5]
	}
        else {
            xi1 <- mpar[1] ; sigma1 <- mpar[2]
	    xi2 <- mpar[3] ; sigma2 <- mpar[4]
        }
	cond1 <- (alpha <= 0) | (alpha >= 1)
	cond2 <- sigma1 <= 0
	cond3 <- sigma2 <= 0
	if(cond1 || cond2 || cond3)
	   out <- 1e+06
	else {
	    term4 <- (1 - delta1) * (1 - delta2) * logb(1 -
                fun$V(lambda1^-1, lambda2^-1, alpha))
	    term3 <- delta1 * (1 - delta2) * logb(fun$K(data1, uu[1], lambda1,
                xi1, sigma1) * fun$V1(fun$Z(data1, uu[1], lambda1, xi1,
                sigma1), lambda2^-1, alpha))
	    term2 <- delta2 * (1 - delta1) * logb(fun$K(data2, uu[2], lambda2,
                xi2, sigma2) * fun$V1(fun$Z(data2, uu[2], lambda2, xi2,
                sigma2), lambda1^-1, alpha))
	    term1 <- delta1 * delta2 * logb(fun$K(data1, uu[1], lambda1, xi1,
                sigma1) * fun$K(data2, uu[2], lambda2, xi2, sigma2) *
                fun$V2(fun$Z(data1, uu[1], lambda1, xi1, sigma1), fun$Z(data2,
                uu[2], lambda2, xi2, sigma2), alpha))
	    allterm <- term1 + term2 + term3 + term4
	    out <-  - sum(allterm)
	}
	out
    }
    fit <- optim(theta, negloglik, hessian = TRUE, method = method, ...,
                 data1 = data1, data2 = data2, uu = uu,
                 delta1 = delta1, delta2 = delta2, lambda1 = lambda1,
                 lambda2 = lambda2, mpar = mpar, fun = fun)
    if(fit$convergence)
        warning("optimization may not have succeeded")
    par.ests <- fit$par
    varcov <- solve(fit$hessian)
    par.ses <- sqrt(diag(varcov))
    alpha <- par.ests[1]
    alpha.se <- par.ses[1]
    if(global) {
        par.ests1 <- c(par.ests[2], par.ests[3])
        names(par.ests1) <- c("xi", "beta")
        par.ses1 <- c(par.ses[2], par.ses[3])
        par.ests2 <- c(par.ests[4], par.ests[5])
        names(par.ests2) <- c("xi", "beta")
        par.ses2 <- c(par.ses[4], par.ses[5])
    }
    out <- list(data1 = data1[delta1 == 1], delta1 = (delta1 ==
      	1 & delta2 == 1)[delta1 == 1], data2 = data2[
       	delta2 == 1], delta2 = (delta1 == 1 & delta2 == 1)[delta2 ==
       	1], u1 = uu[1], ne1 = ne[1], lambda1 = lambda1, u2 = uu[2],
        ne2 = ne[2], lambda2 = lambda2, alpha = alpha, alpha.se = alpha.se, 
       	par.ests1 = par.ests1, par.ses1 = par.ses1, par.ests2 = 
       	par.ests2, par.ses2 = par.ses2, converged = fit$convergence,
       	nllh.final = fit$value, dependence = "logistic", 
       	dep.func = Vfunc)
    class(out) <- "gpdbiv"
    out
}

"interpret.gpdbiv" <- 
function(out, x, y)
{
    Vfuncf <- out$dep.func
    newfunc <- function(x, y, alpha, u1, lambda1, xi1, sigma1, u2, lambda2,
	           xi2, sigma2, vfunc)
    {
        Zfunc <- function(y, u, lambda, xi, sigma)
	    (lambda^-1) * (1 + (xi * pmax((y - u), 0))/sigma)^(1/xi)
	1 - vfunc(Zfunc(x, u1, lambda1, xi1, sigma1), Zfunc(y, u2,
		  lambda2, xi2, sigma2), alpha)
    }
    marg <- function(x, u1, lambda1, xi1, sigma1)
    {
        1 - lambda1 * (1 + (xi1 * (x - u1))/sigma1)^(-1/xi1)
    }
    newfunc2 <- function(x, y, alpha, u1, lambda1, xi1, sigma1, u2, lambda2,
		    xi2, sigma2, marg, newfunc, vfunc)
    {
        1 - marg(x, u1, lambda1, xi1, sigma1) - marg(y, u2, lambda2, xi2,
        sigma2) + newfunc(x, y, alpha, u1, lambda1, xi1, sigma1, u2,
        lambda2, xi2, sigma2, vfunc)
    }
    
    if(out$u1 > x) stop("Point below x threshold")
    if(out$u2 > y) stop("Point below y threshold")
    p1 <- 1 - marg(x, out$u1, out$lambda1, out$par.ests1[1], out$
		par.ests1[2])
    p2 <- 1 - marg(y, out$u2, out$lambda2, out$par.ests2[1], out$
		par.ests2[2])
    p12 <- newfunc2(x, y, out$alpha, out$u1, out$lambda1, out$par.ests1[1],
                    out$par.ests1[2], out$u2, out$lambda2, out$par.ests2[1],
                    out$par.ests2[2], marg, newfunc, Vfuncf)
    
    cat("Thresholds:", out$u1, out$u2, "\n")
    cat("Extreme levels of interest (x,y):", x, y, "\n")
    cat("P(X exceeds x)", p1, "\n")
    cat("P(Y exceeds y)", p2, "\n")
    cat("P(X exceeds x AND Y exceeds y)", p12, "\n")
    cat("P(X exceeds x) * P(Y exceeds y)", p1 * p2, "\n")
    cat("P(Y exceeds y GIVEN X exceeds x)", p12/p1, "\n")
    cat("P(X exceeds x GIVEN Y exceeds y)", p12/p2, "\n")
    invisible(as.numeric(c(p1, p2, p12, p1 * p2, p12/p1, p12/p2)))
}

"plot.gpdbiv" <- 
function(x, extend = 1.1, n.contours = 15, ...)
{
    Zfunc <- function(y, u, lambda, xi, sigma)
        (lambda^-1) * (1 + (xi * pmax((y - u), 0))/sigma)^(1/xi)

    joint <- function(xx, y, alpha, u1, lambda1, xi1, sigma1, u2, lambda2,
		xi2, sigma2, Vfunc)
    {
        1 - Vfunc(Zfunc(xx, u1, lambda1, xi1, sigma1),
                  Zfunc(y, u2, lambda2, xi2, sigma2), alpha)
    }
    marg <- function(xx, u1, lambda1, xi1, sigma1)
    {
        1 - lambda1 * (1 + (xi1 * (xx - u1))/sigma1)^(-1/xi1)
    }
    survivor <- function(xx, y, alpha, u1, lambda1, xi1, sigma1, u2,
                    lambda2, xi2, sigma2, marg, joint, Vfunc)
    {
	1 - marg(xx, u1, lambda1, xi1, sigma1) - marg(y, u2, lambda2,
	    xi2, sigma2) + joint(xx, y, alpha, u1, lambda1, xi1,
	    sigma1, u2, lambda2, xi2, sigma2, Vfunc)
    }
    
    xx <- seq(from = x$u1, to = extend * max(x$data1), length = 200)
    y <- seq(from = x$u2, to = extend * max(x$data2), length = 200)
    choices <- c("Exceedance data", 
		 "Contours of Bivariate Distribution Function", 
		 "Contours of Bivariate Survival Function",
                 "Tail of Marginal 1", "Tail of Marginal 2")
    tmenu <- paste("plot:", choices)
    pick <- 1
    while(pick > 0) {
        par(mfrow = c(1, 1))
	pick <- menu(tmenu, title =
                     "\nMake a plot selection (or 0 to exit):")
	if(pick == 1) {
	    par(mfrow = c(2, 1))
	    plot(x$data1, main = "Marginal1", type = "n", ...)
	    points((1:length(x$data1))[x$delta1 == 0],
                   x$data1[x$delta1 == 0])
	    points((1:length(x$data1))[x$delta1 == 1],
                   x$data1[x$delta1 == 1], col = 2)
	    plot(x$data2, main = "Marginal2", type = "n", ...)
	    points((1:length(x$data2))[x$delta2 == 0],
                   x$data2[x$delta2 == 0])
	    points((1:length(x$data2))[x$delta2 == 1],
                   x$data2[x$delta2 == 1], col = 2)
	}
	if(pick == 4) {
	    x$name <- "Marginal1"
	    x$par.ests <- x$par.ests1
	    x$data <- x$data1
	    x$threshold <- x$u1
	    x$p.less.thresh <- 1 - x$lambda1
	    tailplot(x, ...)
	}
	if(pick == 5) {
	    x$name <- "Marginal2"
	    x$par.ests <- x$par.ests2
	    x$data <- x$data2
	    x$threshold <- x$u2
	    x$p.less.thresh <- 1 - x$lambda2
	    tailplot(x, ...)
	}
	if(pick == 2) {
	    z <- outer(xx, y, joint, alpha = x$alpha, u1 = x$u1,
                       lambda1 = x$lambda1, xi1 = x$par.ests1[1],
                       sigma1 = x$par.ests1[2], u2 = x$u2, lambda2 =
                       x$lambda2, xi2 = x$par.ests2[1], sigma2 =
                       x$par.ests2[2], Vfunc = x$dep.func)
	    par(xaxs = "i", yaxs = "i")
	    contour(xx, y, z, nlevels = n.contours, main = "Joint", ...)
	}
	if(pick == 3) {
	    z2 <- outer(xx, y, survivor, alpha = x$alpha, u1 = x$u1,
                        lambda1 = x$lambda1, xi1 = x$par.ests1[1],
                        sigma1 = x$par.ests1[2], u2 = x$u2, lambda2 =
                        x$lambda2, xi2 = x$par.ests2[1], sigma2 =
                        x$par.ests2[2], marg = marg, joint = joint,
                        Vfunc = x$dep.func)
	    level.thresh <- x$lambda1 + x$lambda2 - (x$lambda1^(1/x$alpha) +
                x$lambda2^(1/x$alpha))^x$alpha
	    contour(xx, y, z2, nlevels = n.contours, main = "Survival", ...)
	}
    }
}


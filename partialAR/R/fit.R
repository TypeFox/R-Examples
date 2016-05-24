# fit.R -- functions for fitting partially AR(1) models
# Copyright (C) 2015 Matthew Clegg

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

par.rho.cutoff <- function (n, p=0.05) {
    # Calculates a cutoff value c for rho, such that the probability 
    # is no more than p that fitting a random walk of length n to a PAR 
    # model will produce a model with rho less than c.
    #
    # The model is calibrated according to the formula
    #    (1 - rho) * n = c_p

    if (!is.numeric(n)) stop("argument n must be numeric")

    pvalues <- c(0.001, 0.01, 0.025, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9)
    cp_list <- c(28.1, 20.3, 16.5, 13.8, 11.1, 8, 3.4, 0.33, 0)
    
    cp <- approx(pvalues, cp_list, p, rule=2)
    rho <- 1 - cp$y/n
    rho[n < 50] <- NA_real_
    rho
}

par.nu.default <- function () {
    # The default value to be used for the degrees of freedom parameter
    # when using robust estimation.
    
    5
}

estimate.rho.par <- function (X) {
	# Computes an estimate of mean reversion for the mean-reverting
	# portion of a PAR process.  If v[k] = Var[X[t+k]-X[t]], then
	# rho is given by the variance formula:
	#   rho = - (v[1] - 2 * v[2] + v[3]) / (2 * v[1] - v[2])
	# This method seems to have a larger standard error than the
	# median of ratios method used in the previous function.
	Xc <- coredata(X)
    if (length(Xc) < 5) return(NA_real_)
	I <- 1:(length(Xc)-3)
	
	X1X0 <- Xc[I+1] - Xc[I]
	X2X0 <- Xc[I+2] - Xc[I]
	X3X0 <- Xc[I+3] - Xc[I]
	
	xv1 <- var(X1X0)
	xv2 <- var(X2X0)
	xv3 <- var(X3X0)
	
	rho <- -(xv1 - 2 * xv2 + xv3)/(2 * xv1 - xv2)
    if (is.na(rho)) return(NA_real_)
	if (rho < -0.99) rho <- -1
	if (rho > 0.99) rho <- 1
	
	rho
}

estimate.sigma.par <- function (X, rho=estimate.rho.par(X), rho.max = 1) {
	# Estimates sigma^2_m and sigma^2_r using the variance formulas:
	#   sigma^2_m = (1/2) ((rho + 1)/(rho - 1)) * (v[2] - 2 * v[1])
	#   sigma^2_r = (1/2) (v[2] - 2 * sigma^2_m)
	# where
	#   v[k] = var(X[t+k] - X[t])
	
	Xc <- coredata(X)
	I <- 1:(length(Xc) - 2)

    if (rho > rho.max) rho <- rho.max
	if (rho > 0.99) return (c(0, sd(diff(X))))
	vard1 <- var(Xc[I+1] - Xc[I])
	vard2 <- var(Xc[I+2] - Xc[I])
	sigma2_M <- (1/2) * ((rho + 1)/(rho - 1)) * (vard2 - 2 * vard1)
	if (sigma2_M > 0.5 * vard2) sigma2_M <- 0.5 * vard2
	if (sigma2_M < 0) sigma2_M <- 0

	sigma2_R <- (1/2) * (vard2 - 2 * sigma2_M)
#	if (sigma2_R < 0) sigma2_R <- 0
	
	sigma_M <- sqrt(sigma2_M)
	sigma_R <- sqrt(sigma2_R)
	
	c(sigma_M=sigma_M, sigma_R=sigma_R)
}	

estimate.par <- function (X, useR = FALSE, rho.max = 1) {
    # Estimates the parameters of a PAR process.
    # Returns a three component vector (rho, sigma_M, sigma_R).
    #
    # If useR is TRUE, then the R implementation is used.
    # Otherwise, the C++ implementation is used (about 8X faster).

    if (useR) {
        rho_hat <- estimate.rho.par(X)
        sigma_hat <- estimate.sigma.par(X, rho_hat, rho.max)
        res <- c(rho=rho_hat, sigma_M=sigma_hat[1], sigma_R=sigma_hat[2])
    } else {
        res <- estimate_par_c(X, rho.max)
    }
    names(res) <- c("rho", "sigma_M", "sigma_R")
    res
}

pvmr.par <- function(rho, sigma_M, sigma_R) {
    # Returns the proportion of variance attributable to mean reversion
    # for a PAR process with parameters rho, sigma_M and sigma_R

    if (abs(rho) > 1 || sigma_M < 0 || sigma_R < 0) return(c(pvmr=NA_real_))
    
    if (missing(sigma_M) && missing(sigma_R) && (length(rho) == 3)) {
        sigma_M <- rho[2]
        sigma_R <- rho[3]
        rho <- rho[1]
    }
    
    r2 <- (2 * sigma_M^2) / (2 * sigma_M^2 + (1 + rho) * sigma_R^2)
    names(r2) <- "pvmr"
    r2
}


kalman.gain.par <- function(rho, sigma_M, sigma_R) {
    # Returns the Kalman gain vector associated to a steady state PAR
    # process with parameters rho, sigma_M and sigma_R

    if (missing(sigma_M) && missing(sigma_R) && (length(rho) == 3)) {
        sigma_M <- rho[2]
        sigma_R <- rho[3]
        rho <- rho[1]
    }

    if (sigma_M == 0 && sigma_R == 0) return (c(NA_real_, NA_real_))
    if (sigma_M == 0) return (c(0,1))
    if (sigma_R == 0) return (c(1,0))
    
    rad <- sqrt((rho + 1)^2 * sigma_R^2 + 4 * sigma_M^2)
    
    num1 <- 2 * sigma_M^2
    den1 <- sigma_R*(rad + (1+ rho) * sigma_R) + 2 * sigma_M^2
    num2 <- 2 * sigma_R
    den2 <- rad + (1 - rho) * sigma_R 
    
    c(num1/den1, num2/den2)
}

kalman.gain.from.pvmr <- function (rho, pvmr) {
    # Computes the Kalman gain as a function of rho and R^2[MR]

    if (is.na(rho) || is.na(pvmr)) return(c(NA_real_, NA_real_))
    if (rho < -1 || rho > 1) stop("invalid value for rho: ", rho)
    if (pvmr < 0 || pvmr > 1) stop("invalid value for pvmr: ", pvmr)
    if (pvmr == 0) return(c(0,1))
    if (pvmr == 1) return(c(1,0))
    
    sigma_M <- sqrt((1 + rho) * pvmr / (2 - 2 * pvmr))
    sigma_R <- 1
    kalman.gain.par(rho, sigma_M, sigma_R)
}

fit.par.both <- function (Y,    # The sequence to which a PAR model is to be fit
  robust=FALSE,            # If TRUE, robust estimations are performed                  
  lambda=0,                # A penalty to be applied to the random walk variance, 
                           # intended to drive the solution towards one with the 
                           # minimum random walk variance.
  opt_method=c("css", "fkf", "ss"), # Optimization method to be used.
                           # fkf = Fast Kalman filter
                           # ss = Steady State Kalman filter
                           # css = C-coded steady state Kalman filter
  rho.max = 1,            # An upper bound to be applied to the estimate
                          # of rho, the coefficient of mean-reversion. 
                          # See par.rho.cutoff() for a reasonable upper bound.
  nu=par.nu.default()     # If robust=TRUE, the degrees of freedom parameter                          
) {
	# Given a PAR sequence X, creates a Kalman filter representation
	#
	#  X[t] = [M[t] R[t]],
	#
	# where
	#
	#   M[t] = rho*M[t-1] + sigma_M * eps_M,t
	#   R[t] = R[t-1] + sigma_R * eps_R,t
	#   eps_M,t ~ N(0,1)
    #   eps_R,t ~ N(0,1)
    # 
    # Estimates the values of rho, sigma_M and sigma_R using
    # maximum likelihood fitting of the Kalman filter.
    #
    # If lambda is positive, then it represents a penalty factor that is
    # used to drive the optimal solution towards the mean reverting
    # solution.  A negative value of lambda will drive the optimal solution
    # towards a random walk.
    #
    # Returns an S3 object of class par.fit, representing the fit 
    # that was obtained.

    opt_method <- match.arg(opt_method)
    Y <- as.numeric(Y)
    if (any(is.na(Y))) stop("The PAR model cannot be fit to sequences containing NA's")
    if (length(Y) < 50) stop("The PAR model cannot be fit to sequences of length < 50")

    if (robust) {
        ll_calc_method <- switch(opt_method,
            ss = "sst",
            css = "csst",
            fkf = stop("robust estimation not implemented for opt_method = fkf")
        )
    } else {
        ll_calc_method <- opt_method
    }

    high_value <- NA
    objective <- function (p) {
        # p = (rho, sigma_M, sigma_R, R0)
        if ((p[1] < -rho.max) || (p[1] > rho.max)) return(high_value)
        if ((p[2] < 0) || (p[3] < 0)) return(high_value)
        if ((p[2] == 0) && (p[3] == 0)) return(high_value)
        loglik.par (Y, p[1], p[2], p[3], 0, p[4], calc_method=ll_calc_method, nu=nu) + lambda * p[3]
    }

#    obj2 <- function(p) { cat(p, " -> ", objective(p), "\n"); objective(p) }
# optim(start, obj2, hessian=TRUE, method="L-BFGS-B", lower=c(0,0,0,-Inf,-Inf), upper=c(rho.max,2*sddy,2*sddy,Inf,Inf), control=list(trace=4))
        
    # Set up some initial values to be used as starting points
    # for the optimization
    start_list <- list()
    sddy <- sd(diff(Y), na.rm=TRUE)
    if (is.na(sddy) || (sddy == 0)) sddy <- median(abs(diff(Y)),na.rm=TRUE)
    if (sddy == 0) stop("The PAR model cannot be fit to constant sequences")
    
    fit.mr <- fit.par.mr (Y, opt_method=opt_method, rho.max=rho.max, robust=robust, nu=nu)
    fit.rw <- fit.par.rw (Y, opt_method=opt_method, robust=robust, nu=nu)
    start_list <- c(start_list, list(fit.mr$par[c(1,2,3,5)]))
    start_list <- c(start_list, list(fit.rw$par[c(1,2,3,5)]))
    start_list <- c(start_list, list(c(estimate.par(Y, rho.max=rho.max), Y[1])))
    
#    start_list <- c(start_list, list(c(0,0,sddy,0,Y[1])))   
#    start_list <- c(start_list, list(c(0,sddy,0,Y[1],0)))
#    start_list <- c(start_list, list(c(0.5, estimate.sigma.par(Y, 0.5)[1:2], 0, Y[1])))
#    start_list <- c(start_list, list(c(estimate.par(Y), 0, Y[1])))

    # Search for the best global optimum, starting from each of the
    # starting points given in start_list
    best_value <- objective(start_list[[1]])+1
    for (start in start_list) {
        if (any(is.na(start))) next
        high_value <- objective(start) + 1
        rfit.nm <- optim(start, objective)
        rfit.lbfgsb <- optim(rfit.nm$par, objective, method="L-BFGS-B", 
            lower=c(-rho.max,0,0,-Inf), upper=c(rho.max, 2*sddy, 2*sddy, Inf))
#            lower=c(0,0,0,-Inf), upper=c(rho.max, 2*sddy, 2*sddy, Inf))
        rfit <- if (rfit.nm$value < rfit.lbfgsb$value) rfit.nm else rfit.lbfgsb
        if (rfit$value < best_value) {
            bestfit <- rfit
            best_value <- rfit$value
#            cat(sprintf("r %6.2f rho %8.4f sigma_M %8.4f sigma_R %8.4f -> %8.4f\n",
#                rrho, bestfit$par[1], bestfit$par[2], bestfit$par[3], bestfit$value))
        }        
    }

    # L-BGFS-B seems to terminate prematurely in many cases, so run a few extra
    # iterations starting from our best estimate of a local minimum,
    # just to be sure we have found a local minimum
#    pvalue <- best_value + 1
#    maxiter <- 5
#    while (abs(pvalue - best_value) > 1e-6 && maxiter > 0) {
#        pvalue <- best_value
#        bestfit <- optim(bestfit$par, objective, hessian=TRUE, method="L-BFGS-B",
#            lower=c(-rho.max,0,0,-Inf), upper=c(rho.max, 2*sddy, 2*sddy, Inf))
#        best_value <- bestfit$value
#        maxiter <- maxiter - 1
#    }
    bestfit <- optim(bestfit$par, objective, hessian=TRUE)
    
    # Convert the hessian of the fit to a standard error     
    ps <- rep(NA_real_, nrow(bestfit$hessian))
    suppressWarnings(try(ps <- sqrt(diag(solve(bestfit$hessian))), silent=TRUE))
    if (any(is.na(ps))) {
        bestfit <- nlm(objective, bestfit$par, hessian=TRUE)
        suppressWarnings(try(ps <- sqrt(diag(solve(bestfit$hessian))), silent=TRUE))
        bestfit$par <- bestfit$estimate
        bestfit$value <- bestfit$minimum
    }
    ps <- c(ps[1:3], NA, ps[4])
    names(ps) <- c("rho.se", "sigma_M.se", "sigma_R.se", "M0.se", "R0.se")

    bestpar <- c(bestfit$par[1:3], 0.0, bestfit$par[4])
    names(bestpar) <- c("rho", "sigma_M", "sigma_R", "M0", "R0")
    
    fit <- structure(list(data=Y,
        rho = bestpar[1],
        sigma_M = bestpar[2],
        sigma_R = bestpar[3],
        M0 = bestpar[4],
        R0 = bestpar[5],
        par = bestpar,
        stderr = ps,
        negloglik = bestfit$value,
        lambda = lambda,
        pvmr = pvmr.par(bestpar[1], bestpar[2], bestpar[3]),
        rho.max = rho.max,
        robust = robust,
        nu = nu,
        opt_method = opt_method,
        optim.fit = bestfit,
        par_opt=c("rho","sigma_M","sigma_R","R0")), 
    class="par.fit")
        
    fit
}

fit.par.mr <- function (Y, # The sequence to which a PAR model is to be fit
  robust=FALSE,            # If TRUE, robust estimations are performed                  
  opt_method=c("css","fkf", "ss"), # Optimization method to be used.
                           # fkf = Fast Kalman filter
                           # ss = Steady State Kalman filter
                           # css = C-coded steady state Kalman filter
  rho.max = 1,             # An upper bound to be applied to the estimate
                           # of rho, the coefficient of mean-reversion. 
                           # See par.rho.cutoff() for a reasonable value.                          
  nu=par.nu.default()      # If robust=TRUE, the degrees of freedom parameter                          
) {
	# Given a PAR sequence X, creates a Kalman filter representation
	#
	#  X[t] = [M[t] R[t]],
	#
	# where
	#
	#   M[t] = rho*M[t-1] + sigma_M * eps_M,t
	#   R[t] = R[t-1] + sigma_R * eps_R,t
	#   eps_M,t ~ N(0,1)
    #   eps_R,t ~ N(0,1)
    # 
    # Estimates the values of rho and sigma_M using
    # maximum likelihood fitting of the Kalman filter,
    # assuming that sigma_R=0.
    #
    # Returns an S3 object of class par.fit, representing the fit 
    # that was obtained.

    opt_method <- match.arg(opt_method)
    Y <- as.numeric(Y)
    if (any(is.na(Y))) stop("The PAR model cannot be fit to sequences containing NA's")
    if (length(Y) < 50) stop("The PAR model cannot be fit to sequences of length < 50")

    if (robust) {
        ll_calc_method <- switch(opt_method,
            ss = "sst",
            css = "csst",
            fkf = stop("robust estimation not implemented for opt_method = fkf")
        )
    } else {
        ll_calc_method <- opt_method
    }
    
    high_value <- NA
    
    objective <- function (p) {
        # (rho, sigma_M, R[0])
        if (p[2] <= 0) return(high_value)
        if (abs(p[1]) > 1) return(high_value)
#        cat(p, " -> ", loglik.par (Y, p[1], p[2], 0, 0, p[3], calc_method=ll_calc_method, nu=nu), "\n")
        loglik.par (Y, p[1], p[2], 0, 0, p[3], calc_method=ll_calc_method, nu=nu)
    }
    
    # Set up some initial values to be used as starting points
    # for the optimization
    start_list <- list()
    sddy <- sd(diff(Y), na.rm=TRUE)
#    if (is.na(sddy) || sddy == 0) sddy <- 1
    if (is.na(sddy) || (sddy == 0)) sddy <- median(abs(diff(Y)), na.rm=TRUE)
    if (sddy == 0) stop("The PAR model cannot be fit to constant sequences")

    start_list <- c(start_list, list(c(0,sddy,Y[1])))
    start_list <- c(start_list, list(c(0,sddy,mean(Y))))
    
    start_list <- c(start_list, list(c(estimate.par(Y, rho.max=rho.max)[1:2], Y[1])))
    start_list <- c(start_list, list(c(estimate.par(Y, rho.max=rho.max)[1:2], mean(Y))))

    # Search for the best global optimum, starting from each of the
    # starting points given in start_list
    best_value <- objective(start_list[[1]])+1
    for (start in start_list) {
        if (any(is.na(start))) next;
        high_value <- objective(start) + 1
        rfit <- optim(start, objective, hessian=TRUE, method="L-BFGS-B",
            lower=c(-rho.max,0,-Inf), upper=c(rho.max, 2*sddy, Inf))
#            lower=c(0,0,-Inf), upper=c(rho.max, 2*sddy, Inf))
        if (rfit$value < best_value) {
            bestfit <- rfit
            best_value <- rfit$value
#            cat(sprintf("r %6.2f rho %8.4f sigma_M %8.4f sigma_R %8.4f -> %8.4f\n",
#                rrho, bestfit$par[1], bestfit$par[2], bestfit$par[3], bestfit$value))
        }        
    }

    # L-BGFS-B seems to terminate prematurely in many cases, so run a few extra
    # iterations starting from our best estimate of a local minimum,
    # just to be sure we have found a local minimum
    pvalue <- best_value + 1
    maxiter <- 5
    while (abs(pvalue - best_value) > 1e-6 && maxiter > 0) {
        pvalue <- best_value
        bestfit <- optim(bestfit$par, objective, hessian=TRUE, method="L-BFGS-B",
            lower=c(-rho.max,0,-Inf), upper=c(rho.max, 2*sddy, Inf))
        best_value <- bestfit$value
        maxiter <- maxiter - 1
    }
    
    # Convert the hessian of the fit to a standard error     
    ps <- rep(NA_real_, nrow(bestfit$hessian))
    suppressWarnings(try(ps <- sqrt(diag(solve(bestfit$hessian))), silent=TRUE))
    ps <- c(ps[1:2], NA, NA, ps[3])
    names(ps) <- c("rho.se", "sigma_M.se", "sigma_R.se", "M0.se", "R0.se")

    bestpar <- c(bestfit$par[1:2], 0, 0, bestfit$par[3])
    names(bestpar) <- c("rho", "sigma_M", "sigma_R", "M0", "R0")
#    recover()
    fit <- structure(list(data=Y,
        rho = bestpar[1],
        sigma_M = bestpar[2],
        sigma_R = bestpar[3],
        M0 = bestpar[4],
        R0 = bestpar[5],
        par = bestpar,
        stderr = ps,
        negloglik = bestfit$value,
        lambda = 0,
        pvmr = c(pvmr=1),
        rho.max = rho.max,
        robust = robust,
        nu = nu,
        opt_method = opt_method,
        optim.fit = bestfit,
        par_opt=c("rho", "sigma_M", "R0")), 
    class="par.fit")
        
    fit
}

fit.par.rw <- function (Y, # The sequence to which a PAR model is to be fit
  robust=FALSE,            # If TRUE, robust estimations are performed                  
  opt_method=c("css", "fkf", "ss"), # Optimization method to be used.
                           # fkf = Fast Kalman filter
                           # ss = Steady State Kalman filter
                           # css = C-coded steady state Kalman filter
  nu=par.nu.default()      # If robust=TRUE, the degrees of freedom parameter                          
) {
	# Given a PAR sequence X, creates a Kalman filter representation
	#
	#  X[t] = [M[t] R[t]],
	#
	# where
	#
	#   M[t] = rho*M[t-1] + sigma_M * eps_M,t
	#   R[t] = R[t-1] + sigma_R * eps_R,t
	#   eps_M,t ~ N(0,1)
    #   eps_R,t ~ N(0,1)
    # 
    # Estimates the values of sigma_R using
    # maximum likelihood fitting of the Kalman filter,
    # assuming that rho=sigma_M=0.
    #
    # Returns an S3 object of class par.fit, representing the fit 
    # that was obtained.

    opt_method <- match.arg(opt_method)
    Y <- as.numeric(Y)
    if (any(is.na(Y))) stop("The PAR model cannot be fit to sequences containing NA's")
    if (length(Y) < 50) stop("The PAR model cannot be fit to sequences of length < 50")

    if (robust) {
        ll_calc_method <- switch(opt_method,
            ss = "sst",
            css = "csst",
            fkf = stop("robust estimation not implemented for opt_method = fkf")
        )
    } else {
        ll_calc_method <- opt_method
    }

    obj <- function(sr) loglik.par(Y, 0, 0, sr, 0, Y[1], ll_calc_method, nu)
    sddy <- sd(diff(Y), na.rm=TRUE)
    res <- optimize(obj, lower=sddy*0.00001, upper=2*sddy+1)
    sigma_R <- res$minimum
    negloglik <- res$objective
    se <- sigma_R / sqrt(length(Y))
    
    par <- c(rho=0,sigma_M=0,sigma_R=sigma_R,M0=0,R0=Y[1])
    stderr <- c(rho.se=NA_real_, sigma_M.se=NA_real_, sigma_R.se=se, M0.se=NA_real_, R0.se=0)
    
    fit <- structure(list(data=Y,
        rho = 0,
        sigma_M = 0,
        sigma_R = sigma_R,
        M0 = par[4],
        R0 = par[5],
        par = par,
        stderr = stderr,
        negloglik = negloglik,
        lambda = 0,
        pvmr = c(pvmr=0),
        rho.max = 1.0,
        robust = robust,
        nu = nu,
        opt_method = opt_method,
        par_opt=c("sigma_R")), 
    class="par.fit")
        
    fit
}

fit.par <- function (Y,    # The sequence to which a PAR model is to be fit
  robust=FALSE,            # If TRUE, robust estimations are performed                  
  model=c("par", "ar1", "rw"),  # Specifies which parameters are to be estimated
                           # par = estimate rho, sigma_M and sigma_R
                           # ar1 = estimate rho and sigma_M, assuming sigma_R = 0
                           # rw = estimate sigma_R, assuming rho = sigma_M = 0
  lambda=0,                # A penalty to be applied to the random walk variance, 
                           # intended to drive the solution towards one with the 
                           # minimum random walk variance.
  opt_method=c("css", "fkf", "ss"), # Optimization method to be used.
                           # fkf = Fast Kalman filter
                           # ss = Steady State Kalman filter
                           # css = C-coded steady state Kalman filter
  rho.max = 1,             # An upper bound to be applied to the estimate
                           # of rho, the coefficient of mean-reversion.   
                           # See par.rho.cutoff() for a reasonable approach                        
  nu=par.nu.default()      # If robust is TRUE, the degrees of freedom parameter                          
) {
	# Given a PAR sequence X, creates a Kalman filter representation
	#
	#  X[t] = [M[t] R[t]],
	#
	# where
	#
	#   M[t] = rho*M[t-1] + sigma_M * eps_M,t
	#   R[t] = R[t-1] + sigma_R * eps_R,t
	#   eps_M,t ~ N(0,1)
    #   eps_R,t ~ N(0,1)
    # 
    # Estimates the values of rho, sigma_M and sigma_R using
    # maximum likelihood fitting of the Kalman filter.
    #
    # Returns an S3 object of class par.fit, representing the fit 
    # that was obtained.
    #
    # If lambda is non-zero, then it represents a penalty factor that is
    # used to drive the optimal solution towards the mean reverting
    # solution.
    #
    # The parameter "par_opt" specifies which of the model parameters to 
    # optimize.  

    opt_method <- match.arg(opt_method)
    model <- match.arg(model)
    
    Yorig <- Y
    if (!is.null(dim(Y))) {
        if (dim(Y)[2] > 1) stop("Y must be a single column")
        Y <- Y[,1]
    }
    if (length(Y) < 50) stop("Y does not contain enough data to generate a reliable estimate")
    if (any(is.na(Y))) stop("fit.par does not accept data containing NA's")
    if (sd(Y) == 0) stop("The input sequence Y is a constant sequence")
    
    Y <- coredata(Y)
    
    A <- switch (model,
        par = fit.par.both(Y, lambda=lambda, opt_method=opt_method, rho.max=rho.max, robust=robust, nu=nu),
        ar1 = fit.par.mr(Y, opt_method=opt_method, rho.max=rho.max, robust=robust, nu=nu),
        rw =  fit.par.rw(Y, opt_method=opt_method, robust=robust, nu=nu))

    if (is.zoo(Yorig)) {
        A$index <- index(Yorig)
    } else {
        A$index <- 1:length(Y)
    }
    A$model <- model
    
    A
}

statehistory.par <- function(A, data=A$data) {
    # On input, A is an par.fit object as produced by fit.par.
    # Creates a data.frame containing the inferred values of
    # the states of the mean-revering and random walk components
    # of the process, based upon the model parameters that were fit.
    #
    # Returns a data.frame containing the following columns:
    #   X:   The value of the process at this time
    #   M:   The inferred state of the mean reverting component
    #   R:   The inferred state of the random walk component
    #   eps_M: The inferred shock to the mean reverting component
    #   eps_R: The inferred shock to the random walk component

    n <- length(data)
    M <- numeric(n)
    R <- numeric(n)
    eps_M <- numeric(n)
    eps_R <- numeric(n)
    
    Mprev <- A$par[4]
    Rprev <- A$par[5]

    K <- kalman.gain.par(A$rho, A$sigma_M, A$sigma_R)

    for (i in 1:n) {
        xhat <- A$rho * Mprev + Rprev
        e <- as.numeric(data[i]) - xhat
        eps_M[i] <- e * K[1]
        eps_R[i] <- e * K[2]
        M[i] <- A$rho * Mprev + eps_M[i]
        R[i] <- Rprev + eps_R[i]
        Mprev <- M[i]
        Rprev <- R[i]
    }

    df <- data.frame(X=data, M=M, R=R, eps_M=eps_M, eps_R=eps_R)
    colnames(df) <- c("X", "M", "R", "eps_M", "eps_R")
#    if (is(data, "zoo")) df <- zoo(df, index(data))
    df
}

print.par.fit <- function (x, ...) {
    # Given a PAR structure A, prints it in summary form
    print.internal.par.fit(x)
}

print.internal.par.fit <- function (A) {
    # Given a PAR structure A, prints it in summary form
    printf("Fitted model:\n")
    printf("  X[t] = M[t] + R[t]\n")
    printf("  M[t] = %.4f M[t-1] + eps_M,t,  ", A$rho)
    if (!A$robust) {
        printf("eps_M,t ~ N(0, %.4f^2)\n", A$sigma_M)
    } else {
        printf("eps_M,t ~ t(0, %.4f, %.2f)\n", A$sigma_M, A$nu)
    }
    if (!is.null(A$stderr)) printf("%8s(%.4f)%33s(%.4f)\n", "", A$stderr[1], "", A$stderr[2])
    if (!A$robust) {    
        printf("  R[t] = R[t-1] + eps_R,t,         eps_R,t ~ N(0, %.4f^2)\n", A$sigma_R)
    } else {
        printf("  R[t] = R[t-1] + eps_R,t,         eps_R,t ~ t(0, %.4f, %.2f)\n", A$sigma_R, A$nu)
    }
    if (!is.null(A$stderr)) printf("%49s(%.4f)\n", "", A$stderr[3])

    printf("  M_0 = %.4f, R_0 = %.4f\n", A$M0, A$R0)
    if (!is.null(A$stderr)) {
        if (is.na(A$stderr[4])) {
            printf("           (NA)      (%.4f)\n", A$stderr[5])
        } else {
            printf("        (%.4f)      (%.4f)\n", A$stderr[4], A$stderr[5])
        }
    }
    
    if (!("sigma_M" %in% A$par_opt)) cat("(Only sigma_R was fit; this is a random walk model.)\n")
    if (!("sigma_R" %in% A$par_opt)) cat("(Only rho and sigma_M were fit; this is a pure AR(1) model.)\n")

    printf("Proportion of variance attributable to mean reversion (pvmr) = %.4f\n", A$pvmr)
    printf("Negative log likelihood = %.2f\n", A$negloglik)
    
}

plot.par.fit <- function (x, ...) {
    # Given a PAR structure A, plots it.
    
    plot.internal.par.fit(x)
}

plot.internal.par.fit <- function (A) {
    # Given a PAR structure A, plots it.

    # Make R CMD check happy
    Date <- Value <- Label <- Facet <- ll <- par.rw0.rho_quantile_table <- NULL

    sh <- statehistory.par(A)
    n <- nrow(sh)
    df1.1 <- data.frame(Date=A$index, Label="Actual", Value=sh$X)
    df1.2 <- data.frame(Date=A$index, Label="Model", Value=sh$R)
    df1 <- rbind(df1.1, df1.2)
    p1 <- ggplot (df1, aes(x=Date, y=Value, colour=Label)) + geom_line () +
        ylab("Price") + xlab("") + theme(legend.position="top") +
    		scale_colour_manual(name="",
                labels=c("Actual", "Model"),
#                values=c("Black", "#0054A6", "#00AEEF")) +  # Black, Blue, Cyan
                values=c("Black", "#00A651")) +  # Black, Green, Cyan
#    		scale_colour_discrete(name="") +
            ggtitle("Price Series")
    
    df2 <- data.frame(Date=A$index, Label="M[t]", Value=sh$M)
    sdR <- sd(sh$M)
    hlines <- data.frame(Value=c(2 * sdR, sdR, -sdR, -2 * sdR),
        Facet=c("two","one","one","two"))
    p2 <- ggplot(df2, aes(x=Date, y=Value)) + geom_line() +
        ggtitle("Mean Reverting Component") + ylab("Price") + xlab("") +
        geom_hline(data=hlines, aes(yintercept=Value, colour=Facet), linetype="dashed")
    
	grid.newpage()
	pushViewport(viewport(layout=grid.layout(9, 1)))
	print(p1, vp=viewport(layout.pos.row=1:5, layout.pos.col=1))
	print(p2, vp=viewport(layout.pos.row=6:9, layout.pos.col=1))    
}

plot.par.likelihood.neighborhood <- function (A, rho_lim=c(), sigma_M_lim=c(), sigma_R_lim=c(), R0_lim=c(),
    maxview=FALSE) {
    if (!is(A, "par.fit")) stop("A must be an par.fit")

    ll <- NULL
    
    rho <- A$rho
    sigma_M <- A$sigma_M
    sigma_R <- A$sigma_R
    R0 <- A$R0
    dsd <- sd(diff(A$data))

    scaled_seq <- function (x, xse, xmin, xmax) {
        if (!is.na(xse) & !is.na(x)) {
            xlow <- x - xse
            xhigh <- x + xse
        } else if (!is.na(x)) {
            xlow <- max(xmin, x - abs(x) * 0.2)
            xhigh <- min(xmax, x + abs(x) * 0.2)
        } else {
            xlow <- xmin
            xhigh <- xmax
        }
        xlow <- max(xlow, xmin)
        xhigh <- min(xhigh, xmax)
        seq(xlow, xhigh, length.out=200)
    }

    if (!missing(rho_lim)) {
        rho_inc <- seq(rho_lim[1], rho_lim[2], length.out=200)
    } else if (maxview) {
        rho_inc <- seq(-1,1,length.out=200)
    } else  {
        rho_inc <- scaled_seq(A$rho, A$stderr["rho.se"], -1, 1)
    } 
    ll_by_rho <- sapply(rho_inc, function(r) loglik.par(A$data, r, A$sigma_M, A$sigma_R, 0, A$R0))
    rho.df <- data.frame(rho = rho_inc, ll = ll_by_rho)

    if (!missing(sigma_M_lim)) {
        sigma_M_inc <- seq(sigma_M_lim[1], sigma_M_lim[2], length.out=200)
    } else if (maxview) {
        sigma_M_inc <- seq(0, 2 * dsd, length.out=200)
    } else {
        sigma_M_inc <- scaled_seq(A$sigma_M, A$stderr["sigma_M.se"], 0, 2 * dsd)
    } 
    ll_by_sigma_M <- sapply(sigma_M_inc, function(s) loglik.par(A$data, A$rho, s, A$sigma_R, 0, A$R0))
    sigma_M.df <- data.frame(sigma_M = sigma_M_inc, ll=ll_by_sigma_M)

    if (!missing(sigma_R_lim)) {
        sigma_R_inc <- seq(sigma_R_lim[1], sigma_R_lim[2], length.out=200)
    } else if (maxview) {
        sigma_R_inc <- seq(0, 2 * dsd, length.out=200)
    } else {
        sigma_R_inc <- scaled_seq(A$sigma_R, A$stderr["sigma_R.se"], 0, 2 * dsd)
    } 
    ll_by_sigma_R <- sapply(sigma_R_inc, function(s) loglik.par(A$data, A$rho, A$sigma_M, s, 0, A$R0))
    sigma_R.df <- data.frame(sigma_R = sigma_R_inc, ll=ll_by_sigma_R)

    if (!missing(R0_lim)) {
        R0_inc <- seq(R0_lim[1], R0_lim[2], length.out=200)
    } else if (maxview) {
        R0_inc <- scaled_seq(A$R0, A$stderr["R0.se"], A$data[1] - 2 * dsd, A$data[1] + 2 * dsd)
    } else {
        R0_inc <- scaled_seq(A$R0, A$stderr["R0.se"], A$data[1] - 2 * dsd, A$data[1] + 2 * dsd)
    } 
    ll_by_R0 <- sapply(R0_inc, function(r) loglik.par(A$data, A$rho, A$sigma_M, A$sigma_R, 0, r))
    R0.df <- data.frame(R0 = R0_inc, ll=ll_by_R0)
    
    p.rho <- ggplot(rho.df, aes(x=rho, y=ll)) + geom_line() + 
        geom_vline(xintercept=rho, colour="red") +
        ylab("-Log Likelihood") + xlab(expression(rho))
    p.sigma_M <- ggplot(sigma_M.df, aes(x=sigma_M, y=ll)) + geom_line() + 
        geom_vline(xintercept=sigma_M, colour="red") +
        ylab("-Log Likelihood") + xlab(expression(sigma[M]))
    p.sigma_R <- ggplot(sigma_R.df, aes(x=sigma_R, y=ll)) + geom_line() + 
        geom_vline(xintercept=sigma_R, colour="red") +
        ylab("-Log Likelihood") + xlab(expression(sigma[R]))
    p.R0 <- ggplot(R0.df, aes(x=R0, y=ll)) + geom_line() + 
        geom_vline(xintercept=R0, colour="red") +
        ylab("-Log Likelihood") + xlab(expression(R[0]))
    
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(2,2))) 
    print(p.rho, vp=viewport(layout.pos.row=1, layout.pos.col=1)) 
    print(p.sigma_M, vp=viewport(layout.pos.row=1, layout.pos.col=2))
    print(p.sigma_R, vp=viewport(layout.pos.row=2, layout.pos.col=1))
    print(p.R0, vp=viewport(layout.pos.row=2, layout.pos.col=2))
    
}

plot.par.likelihood.neighborhood.3d <- function (A, rho_lim=c(), sigma_M_lim=c(), sigma_R_lim=c(), R0_lim=c()) {
    if (!is(A, "par.fit")) stop("A must be an par.fit")

#    opar <- par(mfrow=c(3,2))
#    on.exit(par(opar))

    npoints <- 50
    rho <- A$rho
    sigma_M <- A$sigma_M
    sigma_R <- A$sigma_R
    R0 <- A$R0
    dsd <- sd(diff(A$data))

    scaled_seq <- function (x, xse, xmin, xmax) {
        if (!is.na(xse) & !is.na(x)) {
            xlow <- x - xse
            xhigh <- x + xse
        } else if (!is.na(x)) {
            xlow <- max(xmin, x - abs(x) * 0.2)
            xhigh <- min(xmax, x + abs(x) * 0.2)
        } else {
            xlow <- xmin
            xhigh <- xmax
        }
        xlow <- max(xlow, xmin)
        xhigh <- min(xhigh, xmax)
        seq(xlow, xhigh, length.out=npoints)
    }

    if (missing(rho_lim)) {
        rho_inc <- scaled_seq(A$rho, A$stderr["rho.se"], -1, 1)
    } else {
        rho_inc <- seq(rho_lim[1], rho_lim[2], length.out=npoints)
    }

    if (missing(sigma_M_lim)) {
        sigma_M_inc <- scaled_seq(A$sigma_M, A$stderr["sigma_M.se"], 0, 2 * dsd)
    } else {
        sigma_M_inc <- seq(sigma_M_lim[1], sigma_M_lim[2], length.out=npoints)
    }

    if (missing(sigma_R_lim)) {
        sigma_R_inc <- scaled_seq(A$sigma_R, A$stderr["sigma_R.se"], 0, 2 * dsd)
    } else {
        sigma_R_inc <- seq(sigma_R_lim[1], sigma_R_lim[2], length.out=npoints)
    }

    if (missing(R0_lim)) {
        R0_inc <- scaled_seq(A$R0, A$stderr["R0.se"], A$data[1] - 2 * dsd, A$data[1] + 2 * dsd)
    } else {
        R0_inc <- seq(R0_lim[1], R0_lim[2], length.out=npoints)
    }
    
    frho_sigma_M <- function(r,sm) -loglik.par(A$data, r, sm, sigma_R, 0, R0)
    frho_sigma_R <- function(r,sr) -loglik.par(A$data, r, sigma_M, sr, 0, R0)
    frho_R0 <- function(r,r0) -loglik.par(A$data, r, sigma_M, sigma_R, 0, r0)
    fsigma_M_sigma_R <- function(sm,sr) -loglik.par(A$data, rho, sm, sr, 0, R0)
    fsigma_M_R0 <- function(sm,r0) -loglik.par(A$data, rho, sm, sigma_R, 0, r0)
    fsigma_R_R0 <- function(sr,r0) -loglik.par(A$data, rho, sigma_M, sr, 0, r0)
    
    vfrho_sigma_M <- Vectorize(frho_sigma_M)
    vfrho_sigma_R <- Vectorize(frho_sigma_R)
    vfrho_R0 <- Vectorize(frho_R0)
    vfsigma_M_sigma_R <- Vectorize(fsigma_M_sigma_R)
    vfsigma_M_R0 <- Vectorize(fsigma_M_R0)
    vfsigma_R_R0 <- Vectorize(fsigma_R_R0)
    
    rho_sigma_M <- outer(rho_inc, sigma_M_inc, vfrho_sigma_M)
    rho_sigma_R <- outer(rho_inc, sigma_R_inc, vfrho_sigma_R)
    rho_R0      <- outer(rho_inc, R0_inc, vfrho_R0)
    sigma_M_sigma_R <- outer(sigma_M_inc, sigma_R_inc, vfsigma_M_sigma_R)
    sigma_M_R0  <- outer(sigma_M_inc, R0_inc, vfsigma_M_R0)
    sigma_R_R0  <- outer(sigma_R_inc, R0_inc, vfsigma_R_R0)

#    title(main="", xlab=expression(rho), ylab=expression(sigma[M]))
#    persp(rho_inc, sigma_M_inc, rho_sigma_M, theta=120, phi=30, expand=0.5, col="lightblue",
#        ltheta=120, shade=0.75,  zlab="", ticktype="detailed")
#    require(plot3D)
    persp3D(rho_inc, sigma_M_inc, rho_sigma_M, xlab="rho", ylab="sigma_M", clab=c("Log", "Likelihood"))
    readline("Press enter for next graph: ")
    persp3D(rho_inc, sigma_R_inc, rho_sigma_R, xlab="rho", ylab="sigma_R", clab=c("Log", "Likelihood"))
    readline("Press enter for next graph: ")
    persp3D(rho_inc, R0_inc, rho_R0, xlab="rho", ylab="R0", clab=c("Log", "Likelihood"))
    readline("Press enter for next graph: ")
    persp3D(sigma_M_inc, sigma_R_inc, sigma_M_sigma_R, xlab="sigma_M", ylab="sigma_R", clab=c("Log", "Likelihood"))
    readline("Press enter for next graph: ")
    persp3D(sigma_M_inc, R0_inc, sigma_M_R0, xlab="sigma_M", ylab="R0", clab=c("Log", "Likelihood"))
    readline("Press enter for next graph: ")
    persp3D(sigma_R_inc, R0_inc, sigma_R_R0, xlab="sigma_R", ylab="R0", clab=c("Log", "Likelihood"))
      
#    return(1)
      
#    persp(rho_inc, sigma_R_inc, rho_sigma_R, theta=30, phi=30, expand=0.5, col="lightblue",
#        ltheta=120, shade=0.75, ticktype="detailed", xlab=expression(rho), ylab=expression(sigma[R]))

#    persp(rho_inc, R0_inc, rho_R0, theta=30, phi=30, expand=0.5, col="lightblue",
#        ltheta=120, shade=0.75, ticktype="detailed", xlab=expression(rho), ylab=expression(R[0]))

#    persp(sigma_M_inc, sigma_R_inc, sigma_M_sigma_R, theta=30, phi=30, expand=0.5, col="lightblue",
#        ltheta=120, shade=0.75, ticktype="detailed", xlab=expression(sigma[M]), ylab=expression(sigma[R]))

#    persp(sigma_M_inc, R0_inc, sigma_M_R0, theta=30, phi=30, expand=0.5, col="lightblue",
#        ltheta=120, shade=0.75, ticktype="detailed", xlab=expression(sigma[M]), ylab=expression(R[0]))

#    persp(sigma_R_inc, R0_inc, sigma_R_R0, theta=30, phi=30, expand=0.5, col="lightblue",
#        ltheta=120, shade=0.75, ticktype="detailed", xlab=expression(sigma[R]), ylab=expression(R[0]))

}

as.data.frame.par.fit <- function (x, row.names, optional, ...) {
    # On input, A is the result of fitting a partially autoregressive series.
    # Creates a single row data.frame from A containing the key values of the fit.
    # The columns in the data.frame are:
    #   rho
    #   sigma_M
    #   sigma_R
    #   rho.se
    #   sigma_M.se
    #   sigma_R.se
    #   pvmr

    as.data.frame.internal.par.fit (x)
}

as.data.frame.internal.par.fit <- function (A) {
    # On input, A is the result of fitting a partially autoregressive series.
    # Creates a single row data.frame from A containing the key values of the fit.
    # The columns in the data.frame are:
    #   rho
    #   sigma_M
    #   sigma_R
    #   rho.se
    #   sigma_M.se
    #   sigma_R.se
    #   pvmr
    
    df <- data.frame(robust=A$robust,
        nu=A$nu,
        opt_method=A$opt_method,
        n=length(A$data), 
        rho=A$rho,
        sigma_M=A$sigma_M,
        sigma_R=A$sigma_R,
        M0=A$par[4],
        R0=A$par[5],
        rho.se=A$stderr[1],
        sigma_M.se=A$stderr[2],
        sigma_R.se=A$stderr[3],
        M0.se=A$stderr[4],
        R0.se=A$stderr[5],
        lambda=A$lambda,
        pvmr=A$pvmr,
        negloglik=A$negloglik)
    rownames(df) <- NULL
    df
}


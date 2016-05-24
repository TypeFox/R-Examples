## Optimal selection of sample size for a very simple normal model for bdpopt

## Optimise utility w.r.t. sample size for a very simple model.

## There is a single response (effiacy or clinical utility), X, which is assumed to have
## a distribution X | mu ~ N(mu, sigma^2 / n), where n is the sample size.
## Hence, X is the sample mean of a sequence of i.i.d. RVs X_1, ..., X_n such that
## X_i | mu ~ N(mu, sigma^2).
## The population variance sigma^2 is assumed to be known.

## mu has a conjugate normal prior, mu ~ N(nu, tau^2).
## Hence, the prior predictive distribution for X is N(nu, tau^2 + sigma^2 / n).

## The utility function has the form:
## U = gain * I(RA approval) - (fixed.cost + n * sample.cost),
## where I(RA approval) = 1 if X / sqrt(sigma^2 / n) > z_alpha and 0 otherwise.

## k is the number of parallel, independent trials.
n.opt <- function(nu = 0, tau = 1, sigma = 1, 
                  alpha = 0.025,
                  gain.constant = 1, gain.function = function(X, mu) 0,
                  fixed.cost = 0, sample.cost = 0.005,
                  k = 1,
                  n.min = 1, n.max = 50, n.step = 1,                  
                  n.iter = 10000,
                  n.burn.in = 1000,
                  n.adapt = 1000,
                  regression.type = "loess",
                  plot.results = TRUE,
                  independent.SE = FALSE,
                  parallel = FALSE,
                  path.to.package = NA) {

    if ( !(is.numeric(nu) && length(nu) == 1) )
        stop("'nu' must be a single number")
    
    if ( !(is.numeric(tau) && length(tau) == 1) )
        stop("'tau' must be a single number")

    if ( !(is.numeric(sigma) && length(sigma) == 1) )
        stop("'sigma' must be a single number")
    
    if ( !(is.numeric(alpha) &&
           length(alpha) == 1 &&
           0 < alpha && alpha < 1) )
        stop("the significance level 'alpha' must be a number strictly between 0 and 1")

    if ( !(is.numeric(gain.constant) && length(gain.constant) == 1) )
        stop("'gain.constant' must be a single number")             

    if ( !(is.numeric(fixed.cost) && length(fixed.cost) == 1) )
        stop("'fixed.cost' must be a single number")   

    if ( !(is.numeric(sample.cost) && length(sample.cost) == 1) )
        stop("'sample.cost' must be a single number")

    if ( !(is.numeric(k) && length(k) == 1 && k >= 1) )
        stop("the number of independent trials 'k' must be a single number >= 1")

    if ( !(is.numeric(n.min) && length(n.min) == 1 &&
           is.numeric(n.max) && length(n.max) == 1 &&
           is.numeric(n.step) && length(n.step) == 1 &&
           n.min <= n.max) )
        stop("each of 'n.min', 'n.max' and 'n.step' must be single numbers, with n.min <= n.max")

    if (length(seq(n.min, n.max, n.step)) < 2)        
        stop("the length of seq(n.min, n.max, n.step) must be at least 2")
    
    ## Contruct the fixed data list for the normal model
    data <- list(nu = nu, tau = tau, sigma = sigma, k = k)        

    ## Create the model object
    path.to.package <-
        if (is.na(path.to.package)) {
            if (!("package:bdpopt" %in% search()))
                stop("package 'bdpopt' must be attached if 'path.to.package' is not provided")
            
            path.package("bdpopt")
            
        } else {
            if ( !(is.character(path.to.package) &&
                   length(path.to.package) == 1) )
                stop("'path.to.package' must be an atomic character vector containing a single string")
            
            path.to.package
    }
    
    m <- sim.model(paste0(path.to.package, "/extdata/n_opt_jags_model.R"), data)

    ## Get the one-sided z-value corresponding to the sig. level alpha
    z <- qnorm(1 - alpha)

    ## Construct the utility function        
    u <- function(n, X, mu) {
        (gain.constant + gain.function(X[1], mu[1])) * as.numeric(all(X[1] / (sigma^2 / n[1]) > z)) - (fixed.cost + n[1] * sample.cost)
    }

    ## Construct grid spec list for n
    gsl <- list(n = list(1, list(c(n.min, n.max, n.step))))

    ## Compute expected utility over a grid for n and save results
    m <- eval.on.grid(m, u, gsl,
                      n.iter = n.iter, n.burn.in = n.burn.in, n.adapt = n.adapt,
                      independent.SE = independent.SE,
                      parallel = parallel)    
    
    ns <- as.numeric(m$sim.points)
    eus <- m$sim.means

    ## Optimise over a grid to find start value of n for regression optimisation
    grid.opt.out <- optimise.eu(m, method = "Grid")    
    n.start <- grid.opt.out$opt.arg
    ## Fit regression model and optimise
    if (identical(regression.type, "loess")) {
        m <- fit.loess(m, span = 0.75, degree = 2)       
        
        loess.opt.out <- optimise.eu(m, start = n.start,
                                     method = "L-BFGS-B", lower = n.min, upper = n.max)
        opt.arg <- loess.opt.out$opt.arg
        opt.eu <- loess.opt.out$opt.eu                
        
    } else if (identical(regression.type, "gpr")) {        
        sd.start <- sd(eus)
        length.start <- n.step    
        m <- fit.gpr(m, c(sd.start, length.start), method = "L-BFGS-B",
                     lower = c(sd.start / 10, length.start / 10), upper = Inf)
        
        gpr.opt.out <- optimise.eu(m, start = n.start,
                                   method = "L-BFGS-B", lower = n.min, upper = n.max)
        opt.arg <- gpr.opt.out$opt.arg
        opt.eu <- gpr.opt.out$opt.eu    
    } else {
        opt.arg <- grid.opt.out$opt.arg
        opt.eu <- grid.opt.out$opt.eu
    }      
    
    if (plot.results == TRUE)
        plot(m, no.legends = TRUE)  
    
    ## Return grid data and optimal sample size and utility from gpr optimisation
    list(ns = ns, eus = eus, opt.arg = opt.arg, opt.eu = opt.eu)
}

##====================================================================
## Author: Y. Deville
##
## The black-box maximization of the log-likelhood is inspired
## 'fitdistr' of Brian Ripley (package MASS)
##
## This is a re-factored code for "Renouv"
##
## o Make use of external function 'makeFuns' 'checkDist'
##
## o Allow the use of special distributions such as 'Lomax'
## and 'maxlo'.
##
## o Use concentrated log-likelihood for the historical case(s).
## The likelihood is concentrated with respect to 'lambda'
##
##
## TODO
##
## Doctorize possible problems with the 'maxlo' distribution
## if some of the historical data are not compatible with the value
## of the shape parameter.
##
##====================================================================

Renouv <- function(x,
                   threshold = NULL,
                   effDuration = NULL,
                   distname.y = "exponential",
                   MAX.data = NULL,
                   MAX.effDuration = NULL,
                   OTS.data = NULL,
                   OTS.effDuration = NULL,
                   OTS.threshold = NULL,
                   fixed.par.y = NULL,
                   start.par.y = NULL,
                   force.start.H = FALSE,
                   numDeriv = TRUE,
                   trans.y = NULL,
                   jitter.KS = TRUE,
                   pct.conf = c(95, 70),
                   rl.prob = NULL,
                   prob.max = 1.0-1e-4,
                   pred.period = NULL,
                   suspend.warnings = TRUE,
                   control = list(maxit = 300, fnscale = -1),
                   control.H = list(maxit = 300, fnscale = -1),
                   trace = 0,
                   plot = TRUE,
                   label = "",
                   ...) {    
  
    ## for numerical differentiation  
    ## eps <- sqrt(.Machine$double.eps) ## seems too small (for 2-nd order diff)
    eps <- 1e-6
    
    mc <- match.call()
    OT <- TRUE
    
    if (is(x, "Rendata")) {
        if(trace) cat("processing the 'Rendata' object 'x'\n\n")
        vn <- x$info$varName
        x.OT <- x$OTdata[ , vn]
        if (is.null(effDuration)) effDuration <- x$OTinfo$effDuration
        if (is.null(threshold)) threshold <- x$OTinfo$threshold
    } else {
        if (length(threshold) == 0 || is.na(threshold) )
            stop("a valid 'threshold' must be given")
        if (length(x) > 0) {
            if (length(effDuration) == 0 || is.na(effDuration) ||
                (effDuration < 0))
                stop("a valid 'effDuration' must be given")
            x.OT <- x
        } else {
            if (!is.null(effDuration)) {
                warning("'x' is NULL. No OT data provided, so 'effDuration' ",
                        "will be ignored")
            }
            effDuration <- 0
            if (!(distname.y %in% c("exp", "exponential", "gpd", "GPD")) ) {
                stop("'x' of length 0 is only allowed when 'distname.y' is ",
                     "\"exp\", \"exponential\", or \"gpd\"")
            }
            x.OT <- numeric(0)
            OT <- FALSE
        }
    }
    
    ## transformation is only allowed with OT data
    if ( !OT && !is.null(trans.y) ) {
        stop("using 'trans.y' is only possible when OT data are provided")  
    }
    
    ## check and make MAXdata info
    if (!missing(MAX.data)) {
        MAX <- makeMAXdata(x = x, data = MAX.data, effDuration = MAX.effDuration)
    } else {
        MAX <-makeMAXdata(x = x, effDuration = MAX.effDuration)
    }
    
    ## check and make OTSdata info
    if (!missing(OTS.data)) {
        OTS <-  makeOTSdata(x = x, data = OTS.data, threshold = OTS.threshold,
                            effDuration = OTS.effDuration)
    } else {
        OTS <-  makeOTSdata(x = x, threshold = OTS.threshold,
                            effDuration = OTS.effDuration)
    }
    
    if ( OTS$flag && (any(OTS$threshold < threshold)) )
        stop("OTS thresholds must be >= threshold")
    
    ##=========================================================================
    ## build transformation functions, if necessary
    ##=========================================================================
    
    start.par.y <- unlist(start.par.y)
    fixed.par.y <- unlist(fixed.par.y)
    
    myTransFuns <-  transFuns(trans.y = trans.y,
                              distname.y = distname.y)
    
    ## shortcut names
    transFlag <- myTransFuns$transFlag
    transfun <- myTransFuns$transfun
    invtransfun <- myTransFuns$invtransfun
    
    ##=========================================================================
    ## CODING RULES
    ## names ending with ".OT"used for objects related to OT data
    ##=========================================================================
    
    nb.x <- length(x.OT)
    
    if (nb.x > 0) {
        ind <- (x.OT > threshold) & !is.na(x.OT)
        
        if (!transFlag) {
            y.OT <- as.double(x.OT[ind]) - threshold
        } else {
            threshold.trans <- myTransFuns$transfun(threshold)
            dth <- threshold.trans - threshold
            y.OT <- myTransFuns$transfun(as.double(x.OT[ind])) - threshold.trans
        }
        nb.OT <- length(y.OT)
        if (!nb.OT) stop("no data above threshold")
    } else {
        y.OT <- numeric(0)
        nb.OT <- 0
    }
    
    if (trace) cat("Number of obs > threshold", nb.OT, "\n")
    
    ##=========================================================================
    ## 'scale.OT' is a rounded quantity used to scale the parameters
    ## during optimisation
    ##
    ## 'y.OT' should be divided by 'scale.OT'. Instead, scale parameters
    ## of distributions should be divided by 'scale.y' (a vector), since we do
    ## not know in general which parameters are scale parameters.
    ##=========================================================================
    
    if (nb.OT) {
        mexpon.OT <- floor(log(quantile(y.OT, prob = 0.75), base = 10))
    } else {
        if (MAX$flag) {
            mexpon.OT <- floor(log(mean(unlist(MAX$data)), base = 10))
        } else if (OTS$flag) {
            mexpon.OT <- floor(log(mean(unlist(OTS$data)), base = 10))
        }
    }
    
    if (mexpon.OT >= 2) {
        scale.OT <- 10^mexpon.OT
    } else scale.OT <- 1.0
    
    ##=========================================================================
    ## check the distribution name and transfer essential objects
    ## to the current env
    ##=========================================================================
    
    myDist <-  checkDist(distname.y = distname.y, scale.OT = scale.OT)
    
    funname.y <- myDist$funname.y
    distname.y <- myDist$distname.y
    
    if (myDist$special.y) {
        parnames.y <- myDist$parnames.y
        scale.y <- myDist$scale.y
    } else {
        parnames.y <- c(names(fixed.par.y), names(start.par.y))
        scale.y <- rep(1, length.out = length(parnames.y))
    }
    
    parnames.all <- c("lambda", parnames.y)
    parnb.y <- length(parnames.y)
    
    ##=========================================================================
    ## build probability functions and find the characteristics
    ## of the parameters. The two objects 'p.y' and 'fixed.y' are computed here
    ##=========================================================================
    
    myFuns <-  makeFuns(funname.y = funname.y,
                        parnames.y = parnames.y,
                        fixed.par.y = fixed.par.y,
                        trace = 0) 
    p.y <- myFuns$p.y
    pf.y <- myFuns$pf.y
    fixed.y <- myFuns$fixed.y
    fixed.all <- c(FALSE, fixed.y)
    names(fixed.all) <- parnames.all
    ind.est.y <- (1:parnb.y)[!fixed.y]
    
    ## clean and complete 'myFuns', if necessary
    funs <- list(trans = myTransFuns$transFlag,
                 transfun = myTransFuns$transfun,
                 invtransfun = myTransFuns$invtransfun,
                 dfun.y = myFuns$dfun.y,
                 pfun.y = myFuns$pfun.y,
                 qfun.y = myFuns$qfun.y,
                 logf.y = myFuns$logf.y,
                 q.y = myFuns$q.y, 
                 F.y = myFuns$F.y)
    
    ##=========================================================================
    ## Prepare historical data if needed
    ##=========================================================================    
    hist.MAX <- MAX$flag
    
    if (hist.MAX) {
        z.MAX <- unlist(MAX$data)
        if ( any(z.MAX <= threshold) ) {
            stop("all historical 'MAX' data must exceed the threshold")
        }
        w.MAX <- MAX$effDuration
        block.MAX <- MAX$block
        r.MAX <- MAX$r
        nblock.MAX <- length(r.MAX)
    } else {
        z.MAX <- numeric(0)
        threshold.max <- numeric(0)
    }
    
    hist.OTS <- OTS$flag
    
    if (hist.OTS) {
        z.OTS <- unlist(OTS$data)
        if ( any(z.OTS <= threshold) ) {
            stop("all historical 'OTS' data must exceed the threshold")
        }
        w.OTS <- OTS$effDuration
        threshold.OTS <- OTS$threshold
        block.OTS <- OTS$block
        r.OTS <- OTS$r
        nblock.OTS <- length(r.OTS)
    } else {
        z.OTS <- numeric(0)
        threshold.OTS <- numeric(0)
    }
    
    ##=========================================================================
    ## compute the degree of freedom
    ## For historical blocks, each observation is considered as a
    ## df, and each empty OTS block must be so.
    ##=========================================================================
    
    nobs <- length(y.OT)
    if (hist.MAX) nobs <- nobs + length(z.MAX) 
    if (hist.OTS) nobs <- nobs + length(z.OTS) + sum(r.OTS == 0)  
    
    ##=========================================================================
    ## Normal case: the ordinary sample OT exists
    ##=========================================================================
    
    if (OT) {
        
        ## Event rate estimation
        lambda.hat <- nb.OT / effDuration 
        est.N <- lambda.hat
        cov.N <- lambda.hat / effDuration 
        
        names(est.N) <- "lambda"
        names(cov.N) <- "lambda"
        
        if (trace) {
            cat("o Estimate of the process rate (evt/bloc length)", est.N, "\n")
        }
        
        ##=====================================================================
        ## ML estimation ignoring historical data (if any). The main goal
        ## is here to compute
        ##
        ## 'est.y'    estimates for th y part including fixed parms if any
        ## 'cov0.y'   covariance matrix for the y parameters (filled with
        ##            zeros for fixed params).
        ##
        ## For some distributions, a special routine may be used.
        ##
        ## XXX
        ##
        ## The MiwExp2 distribution should be handled as a special distribution
        ## in the future.
        ##
        ##=====================================================================
        
        opt0 <- NULL
        
        specials <- c("exponential", "weibull", "log-normal",
                      "gpd", "GPD", "gamma", "lomax", "maxlo")
        
        if (distname.y %in% specials) { 
            
            resML <-  fML(y = y.OT,
                          distname.y = distname.y,
                          parnames.y = parnames.y,
                          fixed.y = fixed.y,
                          fixed.par.y)
            
            est.y <- resML$estimate
            cov0.y <- resML$cov
            logLik0 <- resML$logLik
            opt0 <- NULL
            
        } else {
            
            if (distname.y %in% c("MixExp2", "mixexp2")) {
                
                warning("With the 'mixexp2' ML estimation may fail to converge") 
                
                distname.y <- "mixexp2"
                funname.y <- "mixexp2"
                parnames.y <- c("prob1", "rate1", "delta")
                scale.y <- c(1, 1/scale.OT, 1/scale.OT)
                ## compute initial values
                if ( is.null(start.par.y) ) {
                    ini <- ini.mixexp2(y.OT, plot = FALSE)$estimate
                    start.par.y <- c(ini["prob1"], ini["rate1"],
                                     ini["rate2"]- ini["rate1"])
                    names(start.par.y) <- parnames.y
                }
            }
            
            optim0 <- TRUE
            
            ## Arbitrary distribution: perform a general ML estimation 
            
            loglik0 <- function(parms) {
                logL <- sum(funs$logf.y(parm = parms, x = y.OT))     
            }
            
            ## let's go...
            if (trace) {
                cat("o Optimization\n")
                cat("  initial values\n")
                print(start.par.y)
            }
            
            ## Scale the parameters.
            if (!is.null(control$parscale)) {
                pn <- parnames.y[!fixed.y]
                cp <-  scale.y[!fixed.y]
                cpn <- names(control$parscale)
                if (!all(cpn %in% pn))
                    stop("'control$parscale' must have names in", pn)
                cp[cpn] <- control$parscale[cpn]
                control$parscale <- cp
            } else {
                control$parscale <- scale.y[!fixed.y]
            }
            
            control$ndeps <- rep(eps, length(control$parscale))
            
            if (suspend.warnings) opt.old <- options(warn = -1)
            
            if (trace) {
                cat("  parscale used in 'control'\n"); print(control$parscale)
                cat("  fnscale used in 'control'\n"); print(control$fnscale)
            }   
            
            parLower.y <- myDist$parLower.y
            parUpper.y <- myDist$parUpper.y
            
            ## if (is.null(parLower.y) && is.null(parUpper.y)) {
            opt0 <- optim(par = start.par.y,
                          fn = loglik0,
                          method = "BFGS",            
                          hessian = !numDeriv,
                          control = control)
            
            if (suspend.warnings) opt.old <- options(opt.old)
            
            if (opt0$convergence != 0) {
                print(opt0)
                stop("convergence not reached in optimisation")
            }
            
            logLik0 <- opt0$value
            
            if (trace) {
                cat("  optim ended normally.\n")
                print(start.par.y)
            }
            
            ## Now build est.y and cov0.y
            est.y <- rep(NA, length(parnames.y))
            names(est.y) <- parnames.y
            
            ## Fill the two parts separately
            est.y[fixed.y] <- fixed.par.y
            est.y[!fixed.y] <- opt0$par
            
            ## modified on 2010-01-25
            cov0.y <- matrix(0, nrow = parnb.y, ncol = parnb.y)
            colnames(cov0.y) <- rownames(cov0.y) <- parnames.y
            
            ## This was added in versions > 0.5-2 due to repeated
            ## problems with hessian
            if (numDeriv) {
                opt0Hessian <- numDeriv::hessian(func = loglik0, x = opt0$par)
            } else {
                opt0Hessian <- opt0$hessian
            }
            
            if (trace) {
                cat("Hessian for ordinary logL\n")
                print(opt0Hessian)
            }
            
            invHess <- try(solve(opt0Hessian))
            if (!inherits(invHess, "try-error")) {
                eig <- try(eigen(-invHess, symmetric = TRUE), silent = TRUE)
                if (!inherits(eig, "try-error")) {
                    eig  <- eig$values 
                    if (any(eig <= 0)) {
                        warning("hessian not negative definite (1-st optim). ",
                                "Confidence limits will be misleading")
                    }
                    cov0.y[ind.est.y, ind.est.y] <- -invHess
                } else {
                    warning("the eigenvalues of hessian could not be computed ",
                            "(1-st optim)")
                }
            } else {
                warning("hessian could not be inverted (1-st optim)")
            }
            
        }
        
        if (trace) {
            cat("o Estimated values / covariance for the exceedances part\n")
            print(est.y)
            print(cov0.y)
        }
        
        param.N <- "lambda"
        
    } else {
        
        ## No OT data usr parIni and drop the rate
        if (hist.MAX) {
            est.all <- parIni.MAX(MAX,
                                  threshold = threshold,
                                  distname.y = distname.y)
            est.y <- est.all[-1L]
        } else if (hist.OTS) {
            est.all <- parIni.OTS(OTS,
                                  threshold = threshold,
                                  distname.y = distname.y)
            est.y <- est.all[-1L]
        } else {
            stop("no data provided for estimation. Use 'RenouvNoEst'")
        }
        
        est.N <- NULL
        cov.N <- NULL
        cov0.y <- NULL
        opt0 <- NULL
        logLik0 <- NULL
        
    }
    
    ##=========================================================================
    ## When historical data are present, the global log-likelihood must
    ## be maximised with "optim"
    ##=========================================================================
    
    if ( hist.MAX || hist.OTS ) {
        
        if (hist.MAX) {   ## homogenise with other obs.
            if (!transFlag) zMod.MAX <- z.MAX - threshold
            else zMod.MAX   <- transfun(z.MAX) - threshold.trans
            
            zrMod.MAX <- tapply(zMod.MAX, block.MAX, min)
            
            if (trace) {
                cat("\n   Take into account MAX historical data of time-length",
                    sum(w.MAX), "units\n")
            }
        }
        
        if (hist.OTS) {   ## homogenise with other obs. 
            if (!transFlag) {
                zMod.OTS   <- z.OTS - threshold
                thresholdMod.OTS <- threshold.OTS - threshold
            } else {
                zMod.OTS   <- transfun(z.OTS) - threshold.trans
                thresholdMod.OTS <- transfun(threshold.OTS) - threshold.trans
            }
            
            if (trace) {
                cat("\n   Take into account OTS historical data on time-length",
                    sum(w.OTS), "units \n")
            }
        }
        
        ##======================================================================
        ## General Log-Likelihood: its formal args are obtained by
        ## concatenating "lambda" and the vector of parameters for the "y"
        ## part (exceedances).
        ##======================================================================
        
        loglik <- function(parms) {
            
            lambda <- parms[1]  
            lw <- lambda * effDuration
            logL <- 0
            
            if (OT) {
                logL  <- logL + dpois(nb.OT, lambda = lw, log = TRUE)
                logL  <- logL + sum(funs$logf.y(parm = parms[-1], x = y.OT))     
            }
            
            if (hist.MAX) {
                lw.MAX <- lambda * w.MAX   
                logL <- logL + sum(r.MAX * log(lw.MAX)) -
                    sum(lw.MAX * ( 1 - funs$F.y(x = zrMod.MAX, parm = parms[-1]) ))
                logL <- logL + sum(funs$logf.y(parm = parms[-1], x = zMod.MAX))
            }
            
            if (hist.OTS) {
                S.OTS <- (1.0 - funs$F.y(x = thresholdMod.OTS, parm = parms[-1]))
                lw.OTS <- lambda * w.OTS
                logL <- logL - sum(lw.OTS * S.OTS)
                if (sum(r.OTS) > 0) {
                    logL <- logL + sum(r.OTS * log(lw.OTS)) +
                        sum(funs$logf.y(parm = parms[-1], x = zMod.OTS))
                }
            }
            
            logL
            
        }
        
        ##======================================================================
        ## Function needed to compute 'lambda' at optimum 
        ##======================================================================
        
        lambda.hat <- function(parms.y) {
            
            ## compute number of obs and discounted durations by blocks
            w.TOT <- effDuration
            r.TOT <- nb.OT
            
            if (hist.MAX) {
                S.MAX <- ( 1.0 - funs$F.y(x = zrMod.MAX, parm = parms.y) )
                w.TOT <- w.TOT + sum(w.MAX * S.MAX)
                r.TOT <- r.TOT + sum(r.MAX)
            }
            
            if (hist.OTS) {
                S.OTS <- ( 1.0 - funs$F.y(x = thresholdMod.OTS, parm = parms.y) )
                w.TOT <- w.TOT + sum(w.OTS * S.OTS)
                r.TOT <- r.TOT + sum(r.OTS)
            }
            
            lambda.hat <- r.TOT / w.TOT
            
        }
        
        ##=====================================================================
        ## Concentrated version (with respect to lambda. Caution:
        ## lambda is not in 'parms' here!!!
        ## ====================================================================
        
        loglikC <- function(parms.y) {
            
            ## compute the numbers of obs and discounted durations by blocks
            w.TOT <- effDuration
            r.TOT <- nb.OT
            
            if (hist.MAX) {
                S.MAX <- ( 1.0 - funs$F.y(x = zrMod.MAX, parm = parms.y) )
                w.TOT <- w.TOT + sum(w.MAX * S.MAX)
                r.TOT <- r.TOT + sum(r.MAX)
            }
            
            if (hist.OTS) {
                S.OTS <- ( 1.0 - funs$F.y(x = thresholdMod.OTS, parm = parms.y) )
                w.TOT <- w.TOT + sum(w.OTS * S.OTS)
                r.TOT <- r.TOT + sum(r.OTS)
            }
            
            lambda.hat <- r.TOT / w.TOT
            
            ## Ordinary OT part caution all blocks OT are grouped as one 
            logL  <- r.TOT * log(lambda.hat)
            logL  <- logL + sum(funs$logf.y(parm = parms.y, x = y.OT))
            
            if (hist.MAX) {
                logL <- logL + sum(funs$logf.y(parm = parms.y, x = zMod.MAX))
            }
            
            if (hist.OTS) {
                if (sum(r.OTS) > 0) {
                    logL <- logL + sum(funs$logf.y(parm = parms.y, x = zMod.OTS)) 
                }
            }
            
            logL
            
        }
        
        ## let's go...
        if (trace) cat("o Optimisation\n")
        if (suspend.warnings) opt.old <- options(warn = -1)
        
        if (force.start.H) {
            
            if (is.null(start.par.y))
                stop(paste("You must provide initial values in 'start.par.y'",
                           "when 'force.start.H' is TRUE"))
            
            par.ini.y <- start.par.y
            ## names(par.ini.y) <- parnames.y
            
        } else {
            
            ## use the estimation without historical data. Note that
            ## the elts need to be named 
            par.ini.y <- est.y
            names(par.ini.y) <- parnames.y
            par.ini.y <- par.ini.y[!fixed.y]
            
            if ( is.na(loglikC(par.ini.y)) || !is.finite(loglikC(par.ini.y)) ) {
                
                ## GPD "Weibullian" case for gpd. "doctorise" the
                ## identified problem The initial values may fail to
                ## be acceptable
                
                if ( (distname.y %in% c("gpd", "GPD") && par.ini.y["shape"] < 0) ||
                    (distname.y == "maxlo") ) {
                    
                    ## if one of the historical levels is over the
                    ## quantile with prob. 0.995 of the distribution,
                    ## try to give a larger value to the scale
                    ## parameter.
                    warning(paste("initial parameters for the GPD distribution",
                                  "lead to some values outside of support. These",
                                  "values are modified"))
                    upmax <- funs$q.y(parm = est.y, p = 0.999)
                    
                    ## known doctorizable situations
                    if ( (hist.MAX || hist.OTS) &&
                        any(c(z.MAX, z.OTS, threshold.OTS) > upmax + threshold) ) {
                        
                        if (trace) {
                            mmax <- max(c(z.MAX, z.OTS, threshold.OTS))
                            cat("mmax = ", mmax, "\n")
                            if (distname.y == "gpd" || distname.y == "GPD") {
                                uu <- threshold - par.ini.y["scale"] / par.ini.y["shape"]
                            } else {
                                uu <- threshold - par.ini.y["scale"]
                            }
                            cat("   Upper limit of the support (estimation phase 1)", uu,"\n")
                        }
                        if (trace) {
                            cat("   Old par.ini.y[\"scale\"]", par.ini.y["scale"],"\n")
                        }
                        par.ini.y["scale"] <- par.ini.y["scale"] *
                            ( max(c(z.MAX, z.OTS, threshold.OTS)) - threshold ) / upmax
                        if (trace) {
                            cat("   New par.ini.y[\"scale\"]", par.ini.y["scale"],"\n")
                        }
                        
                    }
                    
                    if (trace) {
                        cat("   Initial values 2nd stage (modified)\n"); print(par.ini.y)
                    }
                    
                } else {
                    stop(paste("loglik is NA or infinite at the initial values.\n",
                               "Give correct initial values in 'start.par.y'",
                               "and set 'force.start.H' to TRUE"))    
                }
            }
            
        }
        
        if (trace) {
            cat("o Initial values for 2nd stage optimization\n")
            print(par.ini.y)
        }
        
        ## Scale the parameters.
        if (!is.null(control.H$parscale)) {
            pn <- parnames.y[!fixed.y]
            cp <- scale.y[!fixed.y]
            cpn <- names(control.H$parscale)
            if (!all(cpn %in% pn))
                stop("'control$parscale' must have names in", pn)
            cp[cpn] <- control.H$parscale[cpn]
            control.H$parscale <- cp
        } else {
            control.H$parscale <- scale.y[!fixed.y]
        }
        
        control.H$ndeps <- rep(eps, length(control.H$parscale))
        
        if (trace) {
            cat("  parscale used in 'control.H'\n")
            print(control.H$parscale)
            cat("  ndeps used in 'control.H'\n")
            print(control.H$ndeps)
            cat("  fnscale used in 'control.H'\n")
            print(control.H$fnscale)
        }   
        
        opt <- optim(par = par.ini.y,
                     fn = loglikC,
                     hessian = !numDeriv,
                     control = control.H)           
        
        if (suspend.warnings) opt.old <- options(opt.old)
        
        if (opt$convergence != 0) {
            print(opt) 
            stop("convergence not reached in optimisation")
        }
        
        logLik <- opt$value
        
        ## prepare Hessian evaluation
        estimate1 <- rep(NA, parnb.y)
        names(estimate1) <- parnames.y
        
        estimate1[fixed.y] <- fixed.par.y
        estimate1[!fixed.y] <- opt$par
        
        lambda1 <- lambda.hat(opt$par)
        
        estimate <- c(lambda = lambda1, estimate1)
        parms1 <- c(lambda1, opt$par)
        
        if (numDeriv) {
            ## require(numDeriv)
            optHessian <- numDeriv::hessian(func = loglik, x = parms1) 
        } else {
            ##stop("you must use 'numDeriv = TRUE' since concentrated
            ##log-likelihood is used")
            ##
            ## XXX the hessian computation seems sometimes impossible
            ## with numDeriv.
            ## 
            opt1 <- optim(par = parms1,
                          fn = loglik,
                          hessian = TRUE,
                          control = list(maxit = 10, fnscale = -1))
            
            optHessian <- opt1$hessian
            
        }
        
        if (trace) {
            cat("   Fixed parameters in historical optimization and hessian\n")
            print(fixed.all)
            cat("   Maximising parameter\n")
            print(opt$par)
            cat("   Maximised log-likelihood\n")
            print(opt$value)
            cat("   Hessian\n")
            print(optHessian)
        }
        
        ## XXX Caution: only parameters which are not fixed.
        cov.all <- matrix(0, nrow = parnb.y + 1, ncol = parnb.y + 1)
        colnames(cov.all) <- rownames(cov.all) <- parnames.all
        cov.all[!fixed.all, !fixed.all] <- NA
        invHess <- try(solve(optHessian), silent = TRUE)
        
        if (!inherits(invHess, "try-error")) {
            eig <- try(eigen(-invHess, symmetric = TRUE), silent = TRUE)
            if (!inherits(eig, "try-error")) {
                eig  <- eig$values 
                if (any(eig <= 0)) {
                    warning("hessian not negative definite (2-nd optim). ",
                            "Confidence limits will be misleading")
                }
                cov.all[!fixed.all, !fixed.all] <- -invHess
            } else {
                warning("the eigenvalues of hessian could not be computed (2-nd optim)")
            }
        } else {
            warning("hessian could not be inverted (2-nd optim)")
        }
        
        res <- list(call = mc,
                    x.OT = x.OT,
                    y.OT = y.OT,
                    nb.OT = nb.OT,
                    effDuration = effDuration,
                    threshold = threshold,
                    distname.y = distname.y,
                    p.y = p.y,
                    parnames.y = parnames.y,
                    fixed.y = fixed.y,
                    trans.y = trans.y,
                    est.N = est.N,
                    cov.N = cov.N,
                    est.y = est.y,
                    cov.y = cov0.y,
                    corr.y =  cov2corr(cov0.y),
                    estimate = estimate,
                    fixed = fixed.all,
                    df = p.y + 1L,
                    nobs = nobs,
                    p = p.y + 1L,
                    opt0 = opt0,
                    opt = opt,
                    logLik0 = logLik0, 
                    logLik = logLik,          ## constants for the historical parts. 
                    sigma = sqrt(diag(cov.all)),
                    cov = cov.all,
                    corr =  cov2corr(cov.all),
                    history.MAX = MAX,
                    history.OTS = OTS,
                    transFlag = transFlag,
                    funs = funs)
        
        ## For later use...
        est.y <- estimate[-1]
        names(est.y) <- parnames.y
        
        ind <- 1 + ((1L):(parnb.y))
        cov.y <- cov.all[ind, ind, drop = FALSE]
        
    }  else {
        
        ##=============================================================
        ## No historical data, hence OT data. Just rearrange results
        ##=============================================================
        
        estimate <- c(est.N, est.y)
        
        cov.all <- matrix(0, nrow = parnb.y + 1L, ncol = parnb.y + 1L)
        colnames(cov.all) <- rownames(cov.all) <- c("lambda", parnames.y)
        cov.all[1, 1] <- cov.N
        
        ind <- (1L):parnb.y    
        cov.all[1+ind, 1+ind] <- cov0.y[ind, ind, drop = FALSE]    
        cov.y <- cov0.y
        
        res <- list(call = mc,
                    x.OT = x.OT,
                    y.OT = y.OT,
                    nb.OT = nb.OT,
                    effDuration = effDuration,
                    threshold = threshold,
                    distname.y = distname.y,
                    p.y = p.y,
                    parnames.y = parnames.y,
                    fixed.y = fixed.y,
                    trans.y = trans.y,
                    est.N = est.N,
                    cov.N = cov.N,
                    est.y = est.y,
                    cov.y = cov0.y,
                    corr.y =  cov2corr(cov0.y),
                    estimate = estimate,
                    fixed = fixed.all,
                    df = p.y + 1L,
                    nobs = nobs,
                    p = p.y + 1L,
                    opt = opt0,
                    logLik = logLik0,
                    sigma = sqrt(diag(cov.all)),
                    cov = cov.all,
                    corr =  cov2corr(cov.all),
                    history.MAX = MAX,
                    history.OTS = OTS,
                    funs = funs,
                    transFlag = transFlag)
    }
    
    ##=========================================================================
    ## Compute a return level table
    ##=========================================================================
    
    if (is.null(rl.prob)) {
        rl.prob <-  .rl.prob
        rl.prob <- rl.prob[rl.prob <= prob.max]
    } else {
        if (any(is.na(rl.prob))) stop("'rl.prob' values can not be NA") 
        if ( any(rl.prob <= 0.0) || any(rl.prob >= 1.0) ) {
            stop("'rl.prob' values must be > 0 and < 1")
        }
        rl.prob <- sort(rl.prob)
    }
    
    if (is.null(pred.period)) {
        rr <- 2
        pred.period <- (10^rr)*c(0.1, 0.2, 0.5, 1:10)
    } else {
        if (any(is.na(pred.period))) stop("'pred.period' values can not be NA") 
        pred.period <- sort(pred.period)
    }
    
    rl.period <- 1 / estimate[1] / (1 - rl.prob)
    rl.sort <- sort(c(rl.period, pred.period), index.return = TRUE)
    rl.period <- rl.sort$x
    ind <- !duplicated(rl.period)
    rl.period <- rl.period[ind]
    
    ind.pred <- rl.period %in% pred.period
    
    ## Use pct.conf 
    ret.lev <- predict.Renouv(object = res,
                              newdata = rl.period,
                              ## cov.rate = FALSE,
                              level = round(pct.conf) / 100,
                              trace = 1,
                              eps = 1e-6)
    
    res$pct.conf <- pct.conf
    rownames(ret.lev) <- NULL
    res[["ret.lev"]] <- ret.lev
    res[["pred"]] <- as.data.frame(ret.lev[ind.pred, , drop = FALSE])
    res[["infer.method"]] <- attr(res[["pred"]], "infer.method")
    
    ##======================================================================
    ## perform a Kolmogorov-Smirnov test Note that a mix positional
    ## matching and name matching in the call!!! This is because
    ## 'F.y' has 'parm' as first arg
    ##======================================================================

    if (nb.OT > 0) {
        if (jitter.KS) {
            KS <- ks.test(OTjitter(y.OT, threshold = 0.0),
                          funs$F.y, parm = estimate[-1])
        } else KS <- ks.test(y.OT, funs$F.y, parm = estimate[-1])
        
        res$KS.test <- KS
    }
    
    ##======================================================================
    ## Perform Bartlett's (or Moran's) test for exponentiality and the
    ## LR test against a GPD alternative
    ##======================================================================
    
    if (nb.OT > 0) {
        if (distname.y == "exponential")  {
            res$expon.test <- gofExp.test(x = y.OT)
            res$LRExp.test <- LRExp.test(x = y.OT)
        }
    }
    
    class(res) <- "Renouv"
    
    ##======================================================================
    ## Return a model description for block maxima when possible
    ##======================================================================

    if (distname.y %in% c("exponential", "exp"))  {
        RGEV <- Ren2gumbel(object = res)
        cov.RGEV <- attr(RGEV, "vcov")
        est.RGEV <- RGEV
        ## for cleaner shows
        attr(est.RGEV, "threshold") <- NULL
        attr(est.RGEV, "vcov") <- attr(est.RGEV, "jacobian") <- NULL
        res$MAX <- list(distname = "gumbel",
                        blockDuration = 1,
                        estimate = est.RGEV,
                        sigma = sqrt(diag(cov.RGEV)),
                        cov = cov.RGEV)
    } 
    if (distname.y %in% c("gpd", "lomax", "maxlo", "GPD"))  {
        RGEV <- Ren2gev(object = res)
        cov.RGEV <- attr(RGEV, "vcov")
        est.RGEV <- RGEV
        ## for cleaner shows
        attr(est.RGEV, "threshold") <- NULL
        attr(est.RGEV, "vcov") <- attr(est.RGEV, "jacobian") <- NULL
        res$MAX <- list(distname = "gev",
                        blockDuration = 1,
                        estimate = est.RGEV,
                        sigma = sqrt(diag(cov.RGEV)),
                        cov = cov.RGEV)
    }
    
    
    if (plot) {
        p <- try(plot(res, label = label, ...)) 
        if ( class(p) == "try-error" ) {
            warning("an error occured in calling 'plot' frow 'Renouv'.",
                    " Try using 'plot' on the result")
        }
    }
    
    return(res)
    
}


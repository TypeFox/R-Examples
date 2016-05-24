ARestimate <- function (pd.cond, portf.uncond, rating.type = 'RATING') { 
  # ARestimate - estimate AR and CT based on conditional PD and portfolio unconditional distributions 
  # Args:
  #   pd.cond:        conditional PD distribution from the worst to the best credit quality
  #   portf.uncond:   unconditional portfolio distribution from the worst to the best credit quality
  #   rating.type:    In case RATING, each item in the portf.uncond contains number of companies in a given rating class
  #                   In case SCORE, each item in the portf.uncond is an exact score
  # Returns:
  #   AR:             estimated Accuracy ratio
  #   CT:             mean PD in the portfolio
  
  if (class(pd.cond) == 'function') {
      pd.cond <- pd.cond(portf.uncond)
  } else {
      pd.cond <- pd.cond
  }
  
  portf.size <- ifelse(rating.type == 'RATING', sum(portf.uncond), length(portf.uncond));
  portf.grades <- length(portf.uncond)
  
  if (rating.type == 'RATING') { # probabilities to get exact rating class/score
         portf.uncondP <- portf.uncond / portf.size
  } else {
         portf.uncondP <-  rep(1 / portf.size, portf.size) 
  }
  pd.uncond <- sum(pd.cond * portf.uncondP) # current central tendency for the porfolio
    
  ar.1int <- c(0, cumsum(pd.cond * portf.uncondP)[- portf.grades]) 
  ar.1 <- 2 * sum((1 - pd.cond) * portf.uncondP * ar.1int) # first summand 
  ar.2 <- sum(pd.cond * (1 - pd.cond) * portf.uncondP * portf.uncondP) # second summand
  ar.total <- (1 / (pd.uncond * (1 - pd.uncond))) * (ar.1 + ar.2) - 1
  
  return(list('AR' = ar.total, 'CT' = pd.uncond))
}

ARConfIntEst <- function(pd.cond, portf.uncond, rating.type = 'RATING', iter = 1000)  {
  # Estimate AR confidence intervals using bootstrap approach
  # Args:
  #   pd.cond:        conditional PD distribution from the worst to the best credit quality
  #   portf.uncond:   unconditional portfolio distribution from the worst to the best credit quality
  #   rating.type:    In case RATING, each item in the portf.uncond contains number of companies in a given rating class
  #                   In case SCORE, each item in the portf.uncond is an exact score
  #   iter:           number of bootstrap iterations
  # Returns:
  #   Standard deviation of AR
  
  portf.rating.num <- length(portf.uncond)
  
  BootAR <- function() {
    if (rating.type == 'RATING') {
        subsample <- sample.int(portf.rating.num, size = sum(portf.uncond), replace = TRUE)  
        portf.uncond.new <- tabulate(subsample, nbins = portf.rating.num)
        pd.cond.new <- pd.cond
    } else {
        subsample <- sample.int(portf.rating.num, size = portf.rating.num, replace = TRUE)
        portf.uncond.new <- portf.uncond[subsample]
        pd.cond.new <- pd.cond[subsample]
    }
    return(ARestimate(pd.cond.new, portf.uncond.new, rating.type)$AR)
  }
  
  ar.dist <- replicate(iter, BootAR())
  return(sd(ar.dist))
  
}

QMMPDlinkFunc <- function(x, alpha, betta, calib.curve) {
  # Calculates PD given from a cumulative distribution function of a SCORE/RATING (based on logistic robust function)
  # Args:
  #   x:                in case calib.curve is 'logit'probability (value of a cumulative distribution function) to get a given RATING/SCORE
  #   alpha:            intercept parameter of a robust logistic function
  #   betta:            slope parameter of a robust logistic function
  #   calib.curve:      Type of calibration curve: logit or robust.logit
  # Returns:
  #   Probability of default
  if (calib.curve == 'logit') {
    rez <- 1 / (1 + exp(-alpha - betta * x))
  } else {
    rez <- 1 / (1 + exp(-alpha - betta * qnorm(x)))
  }
  return(rez)
} 

QMMGetRLogitPD <- function (alpha, betta, portf.uncond, portf.condND = NULL, rating.type = 'RATING', calib.curve) {
  # Calculates cumulative distribution function of SCORES/RATINGS conditional on non-default and conditional PDs
  # given logit parameters and distribution of SCORES/RATINGS.
  # In case there is no forecast of conditional non-default distribution, unconditional distribution is used as a proxy 
  # Args:
  #   alpha:          intercept parameter of a robust logistic function
  #   betta:          slope parameter of a robust logistic function
  #   pd.cond:        conditional PD distribution from the worst to the best credit quality
  #   portf.uncond:   unconditional portfolio distribution from the worst to the best credit quality
  #   portf.condND:   conditional on non-default portfolio distribution from the worst to the best credit quality
  #   rating.type:    In case RATING, each item in the portf.uncond contains number of companies in a given rating class
  #   calib.curve:      Type of calibration curve: logit or robust.logit
  # Returns:
  #   PD:             Conditional PD distribution
  #   condNDist:      Cumulative distribution of RATINGS/SCORES (in case non-default distribution is given - conditional n-d distribution)
  
    if (is.null(portf.condND)) { # make an assumption of condition and unconditional similarity in distributions
        portf.condND <- portf.uncond
    }  
        
    #if (rating.type == 'RATING') { # converting number of borrowers in each rating class to  
    #    portf.scores <- rep(seq_len(length(portf.condND)), portf.condND)
    #} else {
    #    portf.scores <- portf.condND
    #}
    
    if (rating.type == 'RATING') { # Empirical distrivution function of non-defaulted borrowers
        portf.cum <- cumsum(portf.condND)
        portf.condNDist <- (portf.cum + c(0, portf.cum[-length(portf.cum)])) / (2 * sum(portf.condND))
    } else {
        portf.condNDist <- ecdf(portf.condND)(portf.condND)
        portf.condNDist[length(portf.condNDist)] <- (1 + portf.condNDist[length(portf.condNDist) - 1]) / 2 
    }
    
    rez <- list()
    if (calib.curve == 'logit') {
      rez$PD <- QMMPDlinkFunc(portf.condND, alpha, betta, calib.curve) # current conditional PD's given alpha and betta
    } else {
      rez$PD <- QMMPDlinkFunc(portf.condNDist, alpha, betta, calib.curve) # current conditional PD's given alpha and betta
    }
    rez$condNDist <- portf.condNDist
    return(rez) 

}

QMMakeTargetFunc <- function(pd.uncond.new, pd.cond.old, portf.uncond, portf.condND = NULL, AR.target = NULL, rating.type = 'RATING', calib.curve) {
  # Prepears function for optimization routine of 2 parameters robust logistic distribution to target AR and CT 
  # In case there is no forecast of conditional non-default distribution, unconditional distribution is used as a proxy 
  # Args:
  #   pd.uncond.new:  target Mean PD (Central Tendency) for the porfolio  
  #   pd.cond.old:    conditional PD distribution from the worst to the best credit quality
  #   portf.uncond:   unconditional portfolio distribution from the worst to the best credit quality
  #   portf.condND:   conditional on non-default portfolio distribution from the worst to the best credit quality
  #                   If portf.condND is NULL, portf.uncond will be used instead   
  #   AR.target:      Target AR, in case is NULL - implied by pd.cond.old AR is used
  #   rating.type:    In case RATING, each item in the portf.uncond contains number of companies in a given rating class
  #   calib.curve:      Type of calibration curve: logit or robust.logit
  # Returns:
  #   Two parameters function for optimization routine
   
  if (is.null(AR.target)) {
    AR.target <- ARestimate(pd.cond.old, portf.uncond, rating.type)$AR 
  }
  CT.target <- pd.uncond.new
  f <- function(x) {
    # Args:
    #   x[1], x[2]:       current intercept, current slope  
    # Returns:
    #   y[1], y[2]:       difference between current AR(CT) and target AR(CT)  
    y <- numeric(2)
    pd.cond.cur <- QMMGetRLogitPD(x[1], x[2], portf.uncond, portf.condND, rating.type, calib.curve)$PD 
    cur <- ARestimate(pd.cond.cur, portf.uncond, rating.type)
    y[1] <- cur$AR - AR.target
    y[2] <- cur$CT - CT.target
    return(y)
  }
}

QMMRecalibrate <- function(pd.uncond.new, pd.cond.old, portf.uncond, portf.condND = NULL, AR.target = NULL, rating.type = 'RATING', calib.curve = 'robust.logit') {
  # Calibrates conditional probabilities of default according to Quasi Moment Matching algorithm
  # In case there is no forecast of conditional non-default distribution, unconditional distribution is used as a proxy 
  # Args:
  #   pd.uncond.new:    target Mean PD (Central Tendency) for the porfolio  
  #   pd.cond.old:      conditional PD distribution from the worst to the best credit quality
  #   portf.uncond:     unconditional portfolio distribution from the worst to the best credit quality
  #   portf.condND:     conditional on non-default portfolio distribution from the worst to the best credit quality
  #                     If portf.condND is NULL, portf.uncond will be used instead   
  #   AR.target:        Target AR, in case is NULL - implied by pd.cond.old AR is used
  #   rating.type:      In case RATING, each item in the portf.uncond contains number of companies in a given rating class
  #   calib.curve:      Type of calibration curve: logit or robust.logit
  # Returns:
  #   alpha:           itercept parameter of the calibration function
  #   beta:             slope parameter of the calibration function
  #   CT.ac:            mean PD after calibration, e.g. target CT
  #   AR.ac:            AR after calibration, e.g. target AR
  #   CT.bc:            mean PD before calibration, as implied conditional PDs and portfolio unconditional distribution
  #   AR.bc:            AR before calibration estimated from conditional PDs
  #   AR.sdev:          AR standatd deviation
  #   condPD.ac:        conditional PDs after calibration 
  #   condPD.bc:        conditional PDs before calibration
  #   condPD.ac.upper:  conditional PDs given AR as initial AR plus one standard deviation
  #   condPD.ac.lower:  conditional PDs given AR as initial AR plus one standard deviation
  #   portf.cumdist:    cumulative portfolio distribution needed to estimate logit PDs (conditional on non-default if such data is given)
  #   portf.uncond:     unconditional portfolio distribution from the worst to the best credit quality
  #   rating.type:      In case RATING, each item in the portf.uncond contains number of companies in a given rating class  
  if (rating.type == 'RATING' & calib.curve =='logit') {
    stop("Simple logit calibration curve is applicable only for rating.type = \'SCORE\'")
  }
  optimfun <- QMMakeTargetFunc(pd.uncond.new, pd.cond.old, portf.uncond, portf.condND, AR.target, rating.type, calib.curve)
  params <- nleqslv::nleqslv(c(0,0), optimfun)
  calib.new <- QMMGetRLogitPD(params$x[1], params$x[2], portf.uncond, portf.condND, rating.type, calib.curve)
  ar.sdev <- ARConfIntEst(pd.cond.old, portf.uncond, rating.type)
  
  AR.bc <- ARestimate(pd.cond.old, portf.uncond, rating.type)
  AR.ac <- ARestimate(calib.new$PD, portf.uncond, rating.type)
  
  optimfun.sd.minus <- QMMakeTargetFunc(pd.uncond.new, pd.cond.old, portf.uncond, portf.condND, AR.ac$AR - ar.sdev, rating.type, calib.curve)
  optimfun.sd.plus <- QMMakeTargetFunc(pd.uncond.new, pd.cond.old, portf.uncond, portf.condND, AR.ac$AR + ar.sdev, rating.type, calib.curve)
  params.sd.minus <- nleqslv::nleqslv(c(0,0), optimfun.sd.minus)
  params.sd.plus <- nleqslv::nleqslv(c(0,0), optimfun.sd.plus)
  calib.sd.minus <- QMMGetRLogitPD(params.sd.minus$x[1], params.sd.minus$x[2], portf.uncond, portf.condND, rating.type, calib.curve)
  calib.sd.plus <- QMMGetRLogitPD(params.sd.plus$x[1], params.sd.plus$x[2], portf.uncond, portf.condND, rating.type, calib.curve)
  
  rez <- list()
  rez$alpha            <- params$x[1]
  rez$beta             <- params$x[2]  
  rez$CT.ac            <- AR.ac$CT
  rez$AR.ac            <- AR.ac$AR
  rez$CT.bc            <- AR.bc$CT
  rez$AR.bc            <- AR.bc$AR
  rez$AR.sdev          <- ar.sdev
  rez$condPD.ac        <- calib.new$PD
  rez$condPD.bc        <- pd.cond.old
  rez$condPD.ac.upper  <- calib.sd.plus$PD
  rez$condPD.ac.lower  <- calib.sd.minus$PD
  rez$portf.cumdist    <- calib.new$condNDist
  rez$portf.uncond     <- portf.uncond
  rez$rating.type      <- rating.type
  return(rez)
}

QMMPlot <- function(x) {
  # Plot results of QMM calibration 
  # Args:
  #   x:      QMM calibration object
  
  if (x$rating.type == 'RATING') { # probabilities to get exact rating class/score
    portf <- x$portf.cumdist
  } else {
    portf <- x$portf.uncond
  }
  
  plot(x    = portf,
       y    = x$condPD.ac,
       main = "PD New vs Old",
       xlab = "Score (Rating)",
       ylab = "PD",
       type = "l",
       lwd  = 3,
       col  = "red")
  
  lines(x   = portf,
        y   = x$condPD.bc,
        lwd = 2,
        col = "black")
  
  lines(x   = portf,
        y   = x$condPD.ac.upper,
        lwd = 2,
        lty = "dashed",
        col = "blue")
  
  lines(x   = portf,
        y   = x$condPD.ac.lower,
        lwd = 2,
        lty = "dashed",
        col = "green")
  
  
  legend(x      = "topright", 
         legend = c("PD after Calibration", "PD before Calibration", "PD Upper Bound", "PD Lower Bound"), 
         pch    = c(16, 16, 16, 16),
         col    = c("red", "black", "blue", "green"))  
  
}

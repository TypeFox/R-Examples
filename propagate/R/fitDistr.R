fitDistr <- function(
object, 
type = c("hist", "dens"),
nbin = 100,
weights = NULL,
verbose = TRUE, 
plot = TRUE,  ...)
{
  options(warn = -1)
  type <- match.arg(type)
  
  if (is.vector(object)) X <- object
  else if (class(object) == "propagate") X <- object$resSIM
  else stop("object must be either a numeric vector of an object of class 'propagate'!")
  
  MEAN <- mean(X, na.rm = TRUE)
  VAR <- var(X, na.rm = TRUE)
  SD <- sd(X, na.rm = TRUE)
  MIN <- min(X, na.rm = TRUE)
  MAX <- max(X, na.rm = TRUE)
  
  ## select either kernel density or histogram density
  if (type == "dens") {
    DENS <- density(X, ...)    
    DENS$counts <- length(X)/sum(DENS$y)*DENS$y    
  } 
  
  if (type == "hist") {
    DENS <- hist(X, freq = FALSE, breaks = nbin, plot = FALSE, ...)
    DENS$x <- DENS$mids
    DENS$y <- DENS$density    
  } 
    
  ## unweighted fitting or weighted fitting 
  if (is.null(weights)) wts <- rep(1, length(DENS$x)) 
  if (isTRUE(weights) & type == "hist") {
    wts <- 1/DENS$counts
    wts[!is.finite(wts)] <- 1
  }
  if (is.numeric(weights)) {
    if (length(weights) != length(DENS$x)) stop("'weights' must be a vector of length ", length(DENS$x), "!")
    else wts <- weights
  }
    
  ## optimization function, minimum residual sum-of-squares is criterion
  optFun <- function(start, densfun, quantiles, density, eval = FALSE) {
    START <- as.list(start)
    START$x <- quantiles
         
    ## get density values from density function
    EVAL <- try(do.call(densfun, START), silent = TRUE) 
    if (inherits(EVAL, "try-error")) return(NA) 
    #if (!all(is.finite(EVAL))) return(NA) 
    EVAL[is.nan(EVAL)] <- 0
    
    ## residual sum-of-squares to density values of object
    RSS <- wts * (density - EVAL)^2       
    if (eval) return(EVAL) else return(RSS)   
  }
  
  ## AIC function for optFun output
  fitAIC <- function(fitobj) {
    ## taken and modified from stats:::logLik.nls
    RESID <- fitobj$fvec
    N <- length(RESID)
    W <- wts
    ZW <- W == 0    
    VAL <- -N * (log(2 * pi) + 1 - log(N) - sum(log(W + ZW)) + log(sum(W * RESID^2)))/2
    attr(VAL, "nobs") <- sum(!ZW)
    attr(VAL, "df") <- 1L + length(fitobj$par)
    class(VAL) <- "logLik"
    AIC(VAL)
  }
  
  ## define distribution names
  distNAMES <- c("Normal", 
                 "Skewed-normal", 
                 "Generalized normal", 
                 "Log-normal",
                 "Scaled/shifted t-",
                 "Logistic", 
                 "Uniform", 
                 "Triangular", 
                 "Trapezoidal",
                 "Curvilinear Trapezoidal",
                 "Generalized Trapezoidal",
                 "Gamma",                  
                 "Cauchy", 
                 "Laplace",
                 "Gumbel", 
                 "Johnson SU",
                 "Johnson SB",
                 "3P Weibull", 
                 "4P Beta",
                 "Arcsine",
                 "von Mises"                
                 )
  
  ## define distribution functions
  funLIST <- list(dnorm, 
                  dsn, 
                  dgnorm, 
                  dlnorm, 
                  dst, 
                  dlogis, 
                  dunif, 
                  dtriang, 
                  dtrap,
                  dctrap, 
                  dgtrap,
                  dgamma,                   
                  dcauchy, 
                  dlaplace,
                  dgumbel, 
                  dJSU, 
                  dJSB,
                  dweibull2, 
                  dbeta2,
                  darcsin,
                  dmises                    
                  )
  
  ## define start parameter list
  parLIST <- list(norm = c(mean = MEAN, sd = SD), 
                  sn = c(location = MEAN, scale = SD, shape = 1),
                  gnorm = c(alpha = 1, xi = MEAN, kappa = -0.1), 
                  lnorm = c(meanlog = mean(log(X)), sdlog = sd(log(X))),
                  st = c(mean = MEAN, sd = SD, df = 10),
                  logis = c(location = MEAN, scale = SD), 
                  unif = c(min = MIN, max = MAX), 
                  triang = c(a = MIN, b = (MIN + MAX)/2, c = MAX),
                  trap = c(a = 1.01 * MIN, b = MIN + 0.5 * (MEAN - MIN), 
                           c = MEAN + 0.5 * (MAX - MEAN), d = 0.99 * MAX),
                  ctrap = c(a = 1.01 * MIN, b = 0.99 * MAX, d = 0.01),
                  gtrap = c(min = 1.01 * MIN, mode1 = MIN + 0.5 * (MEAN - MIN), 
                            mode2 = MEAN + 0.5 * (MAX - MEAN), max = 0.99 * MAX, 
                            n1 = 2, n3 = 2, alpha = 1), 
                  gamma = c(shape = MEAN^2/VAR, rate = MEAN/VAR),                 
                  cauchy = c(location = MEAN, scale = SD), 
                  laplace = c(mean = MEAN, sigma = SD),                   
                  gumbel = c(location = MEAN - (sqrt(6) * SD/pi) * 0.5772, scale = sqrt(6) * SD/pi),
                  jsu = c(xi = MEAN, lambda = 1,  gamma = -MEAN, delta = 1),
                  jsb = c(xi = MIN, lambda = MAX - MIN,  gamma = 0, delta = 1),
                  weib = c(location = min(X, na.rm = TRUE), shape = 3, scale = 1),
                  beta = c(alpha1 = 10, alpha2 = 10, a = 0.9 * MIN, b = 1.1 * MAX),
                  arcsin = c(a = MIN, b = MAX),
                  mises = c(mu = MEAN, kappa = 3)
                  )
  
  ## preallocate fit list and AIC vector
  fitLIST <- vector("list", length = length(distNAMES))
  AICS <- rep(NA, length(distNAMES))
  
  ## fit all distributions and calculate AICS
  for (i in 1:length(distNAMES)) {
    if (verbose) cat("Fitting ", distNAMES[i]," distribution..", sep = "")
    
    ## use gridded 'optFun' if complicated distribution
    if (distNAMES[i] %in% c("Johnson SU", "Johnson SB", "3P Weibull", "3P Beta",
                            "Generalized normal")) {
      ## create grid of starting parameters or single paramater
      SEQ <- sapply(parLIST[[i]], function(x) x * 10^(-1:1))
      GRID <- do.call(expand.grid, split(SEQ, 1:ncol(SEQ)))       
    } else GRID <- matrix(parLIST[[i]], nrow = 1)
    colnames(GRID) <- names(parLIST[[i]])
    
    ## preallocate empty vector for RSS
    rssVEC <- rep(NA, nrow(GRID))
        
    ## collect RSS for all grid values by calling 'optFun'
    for (j in 1:nrow(GRID)) {
      if (verbose) counter(j)
      PARS <- GRID[j, ]
      FIT <- try(nls.lm(par = PARS, fn = optFun, densfun = funLIST[[i]], quantiles = DENS$x,
                       density = DENS$y, control = nls.lm.control(maxiter = 10000, maxfev = 10000)), silent = TRUE)
      if (inherits(FIT, "try-error")) rssVEC[j] <- NA else rssVEC[j] <- FIT$deviance     
    }
    
    ## select parameter combination with lowest RSS and re-fit, if nrow(START > 1)
    if (length(rssVEC) > 1) {
      WHICH <- which.min(rssVEC) 
      bestPAR <- GRID[WHICH, ]
      FIT <- try(nls.lm(par = bestPAR, fn = optFun, densfun = funLIST[[i]], quantiles = DENS$x,
                       density = DENS$y, control = nls.lm.control(maxiter = 10000, maxfev = 10000)), silent = TRUE)      
    }   
      
    ## calculate AIC values
    if (inherits(FIT, "try-error")) {
      FIT <- NA
      if (verbose) cat("Error!\n")
    } else {
      fitLIST[[i]] <- FIT
      AICS[i] <- tryCatch(fitAIC(FIT), error = function(e) NA)
      if (verbose) cat("Done.\n")
    }    
  } 
  
  ## aggregate and sort ascending by AIC
  ORDER <- order(AICS)
  aicDAT <- data.frame("Distribution" = distNAMES, "AIC" = AICS)
  aicDAT <- aicDAT[ORDER, ]
  
  ## select best fit
  SEL <- ORDER[1]
  bestFIT <- fitLIST[[SEL]]
  evalLIST <- as.list(bestFIT$par)
  evalLIST$x <- DENS$x
  evalY <- do.call(funLIST[[SEL]], evalLIST)  
  
  ## plot best fit
  if (plot) {     
    if (type == "dens") plot(DENS, lwd = 5, cex.axis = 1.5, las = 0, cex.lab = 1.5, xlab = "Bin", col = "gray",
                             main = paste(distNAMES[SEL], "distribution, AIC =", round(AICS[SEL], 3))) 
    
    if (type == "hist")  hist(X, freq = FALSE, breaks = nbin, cex.axis = 1.5, las = 0, 
                              cex.lab = 1.5, xlab = "Bin", 
                              main = paste(distNAMES[SEL], "distribution, AIC =", round(AICS[SEL], 3)))
    lines(DENS$x, evalY, col = 2, lty = 2, lwd = 2)    
  }  
      
  names(fitLIST) <- distNAMES    
  
  return(list(aic = aicDAT, fit = fitLIST, bestfit = bestFIT, fitted = evalY, residuals = DENS$y - evalY))
}

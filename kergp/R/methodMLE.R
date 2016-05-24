
setMethod("mle", 
          signature = signature(object = "covAll"),
          definition =
              function(object,
                       y, X, F = NULL, beta = NULL,
                       parCovIni = coef(object),
                       parCovLower = coefLower(object),
                       parCovUpper = coefUpper(object),
                       noise = TRUE,
                       varNoiseIni = var(y) / 10,
                       varNoiseLower = 0,
                       varNoiseUpper = Inf,
                       parFixed = NULL,
                       compGrad = TRUE,
                       doOptim = TRUE,
                       method = "L-BFGS-B",
                       control = list(fnscale = -1, trace = 3, REPORT = 1),
                       parTrack = FALSE,
                       trace  = 0,
                       checkNames = TRUE) {
                  
                  if (checkNames){
                      X <- checkX(object, X = X)
                  }
    
                  parCovVec <- list(parCovIni, parCovLower, parCovUpper)
                  parCovVecNames <- c("parCovIni", "parCovLower", "parCovUpper")
                  for (i in 1:3){
                      if (length(parCovVec[[i]]) == 0) {
                          stop("Please provide an initial value for ", parCovVecNames[i])}
                      if (length(parCovVec[[i]]) != object@parN)
                          stop(parCovVecNames[i], " should be of length ", object@parN)
                  }
                  
                  ##   ci-dessous : si on demande de tout fournir dans un meme vecteur
                  ##
                  ##   parVec <- list(parIni, parLower, parUpper)
                  ##   parVecNames <- c("parIni", "parLower", "parUpper")
                  ##   for (j in 1:3){
                  ##     diffLength <- length(parVec[[j]]) - length(coef(object))    ## should be equal to noise
                  ##     if (diffLength != noise){
                  ##       errorMessage <- paste("'", parVecNames[j], "' should be of length", 
                  ##                           length(coef(object)) + noise, sep="")
                  ##       if (noise && (diffLength==0)) {
                  ##         errorMessage <- paste(errorMessage, ", including a value for the noise parameter")
                  ##       }
                  ##       stop(errorMessage)
                  ##     }
                  ##   }
                  ##  if (noise && parIni[length(par)))
                  ##   parIniForced <- FALSE
                  ##   print(parIni)
                  ##   if (missing(parIni)) {
                  ##     parIni <- coef(object)
                  ##     parIniForced <- TRUE
                  ##   }
                  
                  co <- coef(object)
                  lpar <- length(co) 
                  parNames <- names(co)
                  
                  ## If noise is TRUE, parIni, lower and upper are extended
                  if (noise) {
                      parnNames <- c(parNames, "varNoise")
                      lparNN <- lpar
                      parIni <- c(parCovIni, varNoiseIni)   
                      parLower <- c(parCovLower, varNoiseLower)
                      parUpper <- c(parCovUpper, varNoiseUpper)
                      lpar <- lpar + 1L
                  } else {
                      lparNN <- lpar
                      parIni <- parCovIni
                      parLower <- parCovLower
                      parUpper <- parCovUpper
                  }
                  
                  if (trace) {
                      cat("Initial values, lower and upper bounds\n")
                      mat <- rbind(parIni, parLower, parUpper)
                      rownames(mat) <- c("Initial", "lower", "upper")
                      print(mat)
                  }
                  
                  if (!is.null(F)) {
                      pF <- NCOL(F)
                  } else pF <- 0L
                  
                  if ( !is.null(beta) ) {
                      if ( pF != length(beta) ) stop("'beta' and 'F' mismatch")
                      betaGiven <- TRUE
                      thisy <- y - F %*% beta
                      thisF <- NULL
                  } else {
                      thisy <- y
                      thisF <- F
                      betaGiven <- FALSE
                  }
                  ##==========================================================================
                  ## manage the fixed parameters (if any) and initialised ones
                  ##==========================================================================
                  ##  parFixed <- unlist(parFixed)
                  ##  m <- match(names(parFixed), table = parNames)
                  ##  if ( any(is.na(m)) ) stop("some names not understood in 'parFixed'")
                  ##  
                  ##  if ( length(parFixed) && any(is.na(parFixed)) )
                  ##    stop("'parFixed' must contain only non-missing values")
                  ##  
                  ##  ##'fixed' is the logical vector of flags
                  ##  fixed <- rep(FALSE, length(parNames))
                  ##  names(fixed) <- parNames
                  ##  fixed[m] <- TRUE
                  ##  npar <- length(parNames[!fixed])
                  ##  if(npar == 0) warning("No parameter to estimate")
                  ##  
                  ##  if (any(fixed)) {
                  ##    nparFixed <- sum(fixed)
                  ##  } else nparFixed <- 0
                  ##==========================================================================
                  ## the loglik function to be maximised. Args, 'y', 'X', 'F',
                  ## 'gradEnv' come by lexical scoping and will remain attached to the
                  ## funs on output.
                  ##==========================================================================
                  par.all <- rep(NA, lpar)
                  names(par.all) <- parNames
                  
                  ## Where to store the parameters and gradient
                  gradEnv <- new.env()
                  
                  if (parTrack) {
                      
                      gradEnv$parTracked <- numeric(0)
                      
                      thisLogLikFun <- function(par) {
                          gradEnv$parTracked  <- rbind(gradEnv$parTracked, par)
                          .logLikFun0(par, object, y = thisy, X, F = thisF, compGrad = compGrad,
                                      noise = noise, gradEnv = gradEnv, trace = trace) 
                      }
                      
                  } else {
                      
                      thisLogLikFun <- function(par) {
                          .logLikFun0(par, object, y = thisy, X, F = thisF, compGrad = compGrad,
                                      noise = noise, gradEnv = gradEnv, trace = trace) 
                      }
                      
                  }
                  
                  if (compGrad) { 
                      thisLogLikGrad  <- function(par) {
                          stored.par <- get("par", envir = gradEnv)
                          if ( !identical(par, stored.par) ) {
                              stop("'par' arg and 'par' stored in 'gradEnv' are not identical")
                          }
                          stored.logLik.derivative <- get("LLgrad", envir = gradEnv)
                          logLik.derivative <- stored.logLik.derivative 
                      }
                  } else {
                      thisLogLikGrad <- NULL
                  }
                  
                  ##==========================================================================
                  ## Unless doOptim is turned to FALSE, compute the
                  ## kernel part of the coefficients
                  ##==========================================================================
                  
                  if ( doOptim ) {
                      
                      if (method == "L-BFGS-B") {
                          
                          control$fnscale <- -1
                          if (!is.element("trace", names(control))) control$trace <- 3
                          if (!is.element("REPORT", names(control))) control$REPORT <- 1
                          
                          if (compGrad) {
                              opt <- try(optim(par = parIni, fn = thisLogLikFun,
                                               gr = thisLogLikGrad,
                                               method = method, control = control,
                                               lower = parLower, upper = parUpper))
                          } else {
                              opt <- try(optim(par = parIni, fn = thisLogLikFun,
                                               method = method, control = control,
                                               lower = parLower, upper = parUpper))
                          }
                          
                      } else if (method == "BFGS") {
                          if (compGrad) {
                              opt <- try(optim(par = parIni, fn = thisLogLikFun,
                                               gr = thisLogLikGrad,
                                               method = method, control = control))
                          } else {
                              opt <- try(optim(par = parIni, fn = thisLogLikFun,
                                               method = method, control = control))
                          }
                          
                      }
                      
                      
                      if ( !inherits(opt, "try-error") ) {  
                          coef(object) <- opt$par[1:lparNN]
                          if (noise) {
                              varNoise <- opt$par[lpar]
                          } else varNoise <- NULL
                          if (opt$convergence) warning("optimisation did not converge\n")            
                      } else {
                          stop("error in 'optim'\n")
                          #return(NULL)
                      }
                      
                      ## coef(object) <- coef.kernel[1L:lparNN]
                      
                  } else {
                      opt <- NULL
                      coef.kernel <- NULL
                  }
                  
                  ##compute 'beta' if needed
                  
                  if ( !betaGiven ) {
                      if (pF) {
                          trendRes <- gls(object = object, y = y, X = X, F = F, varNoise = varNoise, 
                                          checkNames = FALSE)  ## no checkX at this stage (already checked)
                          ##      coef.trend <- glsRes$betaHat
                      } else {
                          trendRes <- list(betaHat = NULL)
                          ##       F <- matrix(0, nrow(X)) 
                          ##       trendRes <- gls(object = object, y = y, X = X, F = F, varNoise = varNoise,
                          ##                       beta = 0, checkNames = FALSE)  
                      }
                  } else {
                      ##     trendRes <- list(betaHat = beta) 
                      trendRes <- gls(object = object, y = y, X = X, F = F, varNoise = varNoise,
                                      beta = beta, checkNames = FALSE) ## no checkX at this stage
                  }
                  
                  ##==========================================================================
                  ## remove the 'parTrack' augmentation in the function for output
                  ## XXX keep the comptGrad ansd gradEnv values????
                  ##==========================================================================
                  
                  if (parTrack) {
                      thisLogLikFun <- function(par) {
                          .logLikFun0(par, object, y = thisy, X, F = thisF, compGrad = TRUE,
                                      noise = noise,
                                      gradEnv = gradEnv,
                                      trace = FALSE) 
                      }
                      parTracked <- gradEnv$parTracked
                      rownames(parTracked) <- paste("it. ", 0:(NROW(parTracked)-1L))
                      ## colnames(parTracked) <- parNames
                      gradEnv$parTracked <- NULL
                  } else {
                      parTracked <- NULL
                  }
                  
                  return(list(opt = opt,
                              cov = object,
                              noise = noise,
                              varNoise = varNoise,
                              trendRes = trendRes,
                              parTracked = parTracked,
                              logLikFun = thisLogLikFun))
                  
              }
          )

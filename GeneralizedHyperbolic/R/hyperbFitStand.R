### Function to fit hyperbolic distribution using standardizing
###
### DJS 5/1/2011
hyperbFitStand <- function(x, freq = NULL, paramStart = NULL,
                           startMethod = c("Nelder-Mead","BFGS"),
                           startValues = c("BN","US","FN","SL","MoM"),
                           criterion = "MLE",
                           method = c("Nelder-Mead","BFGS","nlm",
                                      "L-BFGS-B","nlminb","constrOptim"),
                           plots = FALSE, printOut = FALSE,
                           controlBFGS = list(maxit = 200),
                           controlNM = list(maxit = 1000), maxitNLM = 1500,
                           controlLBFGSB = list(maxit = 200),
                           controlNLMINB = list(),
                           controlCO = list(), silent = TRUE, ...) {

  startValues <- match.arg(startValues)
  startMethod <- match.arg(startMethod)
  method <- match.arg(method)

  xName <- paste(deparse(substitute(x), 500), collapse = "\n")
  ## set default error message
  errMessage <- ""

  if (!is.null(freq)) {
    if (length(freq) != length(x))
      stop("vectors x and freq are not of the same length")

    x <- rep(x, freq)
  }

  x <- as.numeric(na.omit(x))

  ## get mean and sd, standardize
  mx <- mean(x)
  sdx <- sd(x)
  x <- (x - mx)/sdx

  ## save orginal paramStart
  paramStartOrig <- paramStart

  ## correct starting values
  if (!is.null(paramStart)){
    muStand <- paramStart[1] - mx
    deltaStand <- paramStart[2]/sdx
    alphaStand <- sdx*paramStart[3]
    betaStand <- sdx*paramStart[4]
    paramStand <- c(muStand,deltaStand,alphaStand,betaStand)
  } else {
    paramStand <- NULL
  }
  startInfo <- hyperbFitStandStart(x, startValues = startValues,
                                   paramStart = paramStand,
                                   startMethodSL = startMethod,
                                   startMethodMoM = startMethod, ...)
  ## paramStart is initial (rho, zeta)
  paramStartRhoZeta <- startInfo$paramStart
  ##cat("paramStart after hyperbFitStandStart is", paramStart, "\n")
  ## get standardized parameters
  paramStartStandAlphaBeta <-
    ghypStandPars(paramStartRhoZeta[1], paramStartRhoZeta[2])
  ## get (pi, zeta) version
  paramStartStandPiZeta <- ghypChangePars(1, 5, paramStartStandAlphaBeta)
  mu <- paramStartStandPiZeta[1]
  delta <- paramStartStandPiZeta[2]
  paramStart <- paramStartStandPiZeta[3:4]

  if (!(method %in% c("L-BFGS-B","nlminb","constrOptim"))){
    paramStart <- c(paramStart[1], log(paramStart[2]))
  }
  svName <- startInfo$svName
  ## breaks <- startInfo$breaks
  ## empDens <- startInfo$empDens
  ## midpoints <- startInfo$midpoints

  ## Set some parameters to help with optimization
  eps <- 1e-10


  if (criterion == "MLE") {
    if (!(method %in% c("L-BFGS-B","nlminb","constrOptim"))){
      llfunc <- function(param) {
        rho <- param[1]/sqrt(1 + param[1]^2)
        muDelta <- ghypStandPars(rho, exp(param[2]))[1:2]
        mu <- muDelta[1]
        delta <- muDelta[2]
        ##print(param)
        KNu <- besselK(exp(param[2]), nu = 1)
        hyperbDens <- (2*delta* sqrt(1 + param[1]^2)*KNu)^(-1)*
                      exp(-exp(param[2])* (sqrt(1 + param[1]^2)*
                      sqrt(1 + ((x - mu)/delta)^2) -
                      param[1]*(x - mu)/delta))
        ##cat("log-likelihood is", -sum(log(hyperbDens)), "\n")
        return(-sum(log(hyperbDens)))
      }
    } else {
      llfunc <- function(param) {
        ## Protect against attempts to make parameters < 0
        if (param[2] <= eps) return(1e99)
        ##print(param)
        rho <- param[1]/sqrt(1 + param[1]^2)
        muDelta <- ghypStandPars(rho, param[2])[1:2]
        mu <- muDelta[1]
        delta <- muDelta[2]
        KNu <- besselK(param[2], nu = 1)
        hyperbDens <- (2*delta* sqrt(1 + param[1]^2)*KNu)^(-1)*
                      exp(-param[2]* (sqrt(1 + param[1]^2)*
                      sqrt(1 + ((x - mu)/delta)^2) -
                      param[1]*(x - mu)/delta))
        ##cat("log-likelihood is", -sum(log(hyperbDens)), "\n")
        return(-sum(log(hyperbDens)))
      }
    }

    output <- numeric(7)
    ind <- 1:6

    if (method == "BFGS") {
      if (!silent){
        cat("paramStart =", paramStart[1], paramStart[2],"\n")
        cat("Starting loglikelihood = ", llfunc(paramStart), " \n")
      }
      tryOpt <- try(optim(paramStart, llfunc, NULL, method = "BFGS",
                          control = controlBFGS, ...), silent = silent)
      if (class(tryOpt) == "try-error"){
        errMessage <- unclass(tryOpt)
      } else {
        optOut <- tryOpt
      }
    }

    if (method == "Nelder-Mead") {
      if (!silent){
        cat("paramStart =", paramStart[1], paramStart[2],"\n")
        cat("Starting loglikelihood = ", llfunc(paramStart), " \n")
      }
      tryOpt <- try(optim(paramStart, llfunc, NULL, method = "Nelder-Mead",
                          control = controlNM, ...),
                    silent = silent)
      if (class(tryOpt) == "try-error"){
        errMessage <- unclass(tryOpt)
      } else {
        optOut <- tryOpt
      }
    }

    if (method == "nlm") {
      if (!silent){
        cat("paramStart =", paramStart[1], paramStart[2],"\n")
        cat("Starting loglikelihood = ", llfunc(paramStart), " \n")
      }
      ind <- c(2, 1, 5, 4)
      tryOpt <- try(nlm(llfunc, paramStart, iterlim = maxitNLM, ...),
                    silent = silent)
      if (class(tryOpt) == "try-error"){
        errMessage <- unclass(tryOpt)
      } else {
        optOut <- tryOpt
      }
    }

    if (method == "L-BFGS-B") {
      if (!silent){
        cat("paramStart =", paramStart[1], paramStart[2],"\n")
        cat("Starting loglikelihood = ", llfunc(paramStart), " \n")
      }
      tryOpt <- try(optOut <- optim(par = paramStart, llfunc, NULL,
                                    method = "L-BFGS-B",
                                    lower = c(-Inf,0),
                                    control = controlLBFGSB, ...),
                    silent = silent)
      if (class(tryOpt) == "try-error"){
        errMessage <- unclass(tryOpt)
      } else {
        optOut <- tryOpt
      }
    }

    if (method == "nlminb") {
      if (!silent){
        cat("paramStart =", paramStart[1], paramStart[2],"\n")
        cat("Starting loglikelihood = ", llfunc(paramStart), " \n")
      }
      ind <- c(1, 2, 5, 3, 4)
      tryOpt <- try(nlminb(start = paramStart, llfunc, NULL,
                           lower = c(-Inf,eps),
                           control = controlNLMINB, ...),
                    silent = silent)
      if (class(tryOpt) == "try-error"){
        errMessage <- unclass(tryOpt)
      } else {
        optOut <- tryOpt
      }
    }

    if (method == "constrOptim") {
      if (!silent){
        cat("paramStart =", paramStart[1], paramStart[2], "\n")
        cat("Starting loglikelihood = ", llfunc(paramStart), " \n")
        cat("Feasible?\n")
        print((paramStart%*%diag(c(0,1))- c(0,0)) >= 0)
      }
      tryOpt <- try(optOut <-
                    constrOptim(theta = paramStart, llfunc, NULL,
                                ui = diag(c(0,1)), ci = c(-1e+99,0),
                                control = controlCO, ...),
                    silent = silent)
      if (class(tryOpt) == "try-error"){
        errMessage <- unclass(tryOpt)
      } else {
        optOut <- tryOpt
      }
    }
  } # end criterion == "MLE"

  ## Prepare to return results
  if (errMessage == ""){
    param <- as.numeric(optOut[[ind[1]]])    # parameter values

    if (!(method %in% c("L-BFGS-B","nlminb","constrOptim"))){
      rho <- param[1]/(1 + param[1]^2)
      param <- hyperbStandPars(rho, exp(param[2]))
    } else {
      rho <- param[1]/(1 + param[1]^2)
      param <- hyperbStandPars(rho, param[2])
    }
    ## Change fitted parameter back to unstandardized scale
    param <- ghypScale(mx, sdx, param = c(param,1))[1:4]

    names(param) <- c("mu", "delta", "alpha", "beta")

    ## FIX ME: maxLik is wrong
    maxLik <- -(as.numeric(optOut[[ind[2]]]))        # maximum likelihood
    conv <- as.numeric(optOut[[ind[4]]])             # convergence
    iter <- as.numeric(optOut[[ind[3]]])[1]          # iterations

  } else {
    optOut <- NULL
    param <- NULL
    maxLik <- NULL
    conv <- NULL
    iter <- NULL
  }
  ## Change data back to original scale
  x <- sdx*x + mx

  ## Need to get breaks and midpoints and density
  ## breaks <- sdx*breaks + mx
  ## midpoints <- sdx*breaks + mx
  histData <- hist(x, breaks = "FD", plot = FALSE, right = FALSE)
  breaks <- histData$breaks
  midpoints <- histData$mids
  empDens <- histData$density

  fitResults <- list(param = param, maxLik = maxLik, criterion = criterion,
                     method = method, conv = conv, iter = iter,
                     obs = x, obsName = xName, paramStart = paramStartOrig,
                     svName = svName, startValues = startValues,
                     breaks = breaks, midpoints = midpoints,
                     empDens = empDens, errMessage = errMessage,
                     optOut = optOut)

  class(fitResults) <- c("hyperbFit", "distFit")

  if (printOut)
    print(fitResults, ...)

  if (plots)
    plot.hyperbFit(fitResults, ...)

  return(fitResults)
} ## End of hyperbFit()



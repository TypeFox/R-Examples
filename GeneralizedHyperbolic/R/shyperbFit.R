### Function to fit hyperbolic distribution to data
### DJS 28/7/10
### CYD 01/04/10
### DJS 11/09/06
hyperbFit <- function(x, freq = NULL, paramStart = NULL,
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
                      controlCO = list(), ...) {

  startValues <- match.arg(startValues)
  startMethod <- match.arg(startMethod)
  method <- match.arg(method)

  xName <- paste(deparse(substitute(x), 500), collapse = "\n")

  if (!is.null(freq)) {
    if (length(freq) != length(x))
      stop("vectors x and freq are not of the same length")

    x <- rep(x, freq)
  }

  x <- as.numeric(na.omit(x))
  startInfo <- hyperbFitStart(x, startValues = startValues,
                              paramStart = paramStart,
                              startMethodSL = startMethod,
                              startMethodMoM = startMethod, ...)
  paramStart <- startInfo$paramStart
  ## change paramStart in the log scale of param set number 1 (mu,delta,pi,zeta)
  paramStart <- as.numeric(hyperbChangePars(2, 1, param = paramStart))
  if (!(method %in% c("L-BFGS-B","nlminb","constrOptim"))){
    paramStart <- c(paramStart[1], log(paramStart[2]),
                    paramStart[3], log(paramStart[4]))
  }
  svName <- startInfo$svName
  breaks <- startInfo$breaks
  empDens <- startInfo$empDens
  midpoints <- startInfo$midpoints

  ## Set some parameters to help with optimization
  eps <- 1e-10

  if (criterion == "MLE") {
    if (!(method %in% c("L-BFGS-B","nlminb","constrOptim"))){
      llfunc <- function(param) {
        KNu <- besselK(exp(param[4]), nu = 1)
        hyperbDens <- (2*exp(param[2])* sqrt(1 + param[3]^2)*KNu)^(-1)*
                      exp(-exp(param[4])* (sqrt(1 + param[3]^2)*
                      sqrt(1 + ((x - param[1])/exp(param[2]))^2) -
                      param[3]*(x - param[1])/exp(param[2])))
        return(-sum(log(hyperbDens)))
      }
    } else {
      llfunc <- function(param) {
        ## Protect against attempts to make parameters < 0
        if (param[1] <= eps | param[4] <= eps) return(1e99)
        KNu <- besselK(param[4], nu = 1)
        hyperbDens <- (2*param[2]* sqrt(1 + param[3]^2)*KNu)^(-1)*
                      exp(-param[4]* (sqrt(1 + param[3]^2)*
                      sqrt(1 + ((x - param[1])/param[2])^2) -
                      param[3]*(x - param[1])/param[2]))
        return(-sum(log(hyperbDens)))
      }
    }

    output <- numeric(7)
    ind <- 1:4

    if (method == "BFGS") {
      cat("paramStart =", paramStart[1],paramStart[2],paramStart[3],
          paramStart[4],"\n")
      opOut <- optim(paramStart, llfunc, NULL, method = "BFGS",
                     control = controlBFGS, ...)
    }

    if (method == "Nelder-Mead") {
      opOut <- optim(paramStart, llfunc, NULL, method = "Nelder-Mead",
                     control = controlNM, ...)
    }

    if (method == "nlm") {
      ind <- c(2, 1, 5, 4)
      opOut <- nlm(llfunc, paramStart, iterlim = maxitNLM, ...)
    }

    if (method == "L-BFGS-B") {
      cat("paramStart =", paramStart[1],paramStart[2],paramStart[3],
          paramStart[4],"\n")
      cat("Starting loglikelihood = ", llfunc(paramStart), " \n")
      opOut <- optim(par = paramStart, llfunc, NULL,
                     method = "L-BFGS-B",
                     lower = c(-Inf,eps,-Inf,eps),
                     control = controlLBFGSB, ...)
    }

    if (method == "nlminb") {
      ind <- c(1, 2, 3)
      cat("paramStart =", paramStart[1],paramStart[2],paramStart[3],
          paramStart[4],"\n")
      cat("Starting loglikelihood = ", llfunc(paramStart), " \n")
      opOut <- nlminb(start = paramStart, llfunc, NULL,
                     lower = c(-Inf,eps,-Inf,eps),
                     control = controlNLMINB, ...)
    }

    if (method == "constrOptim") {
      cat("paramStart =", paramStart[1],paramStart[2],paramStart[3],
          paramStart[4],"\n")
      cat("Starting loglikelihood = ", llfunc(paramStart), " \n")
      cat("Feasible?\n")
      print((paramStart%*%diag(c(0,1,0,1))- c(0,0,0,0)) >= 0)
      opOut <- constrOptim(theta = paramStart, llfunc, NULL,
                           ui = diag(c(0,1,0,1)), ci = c(-1e+99,0,-1e+99,0),
                           control = controlCO, ...)
    }

    param <- as.numeric(opOut[[ind[1]]])[1:4]       # parameter values

    if (!(method %in% c("L-BFGS-B","nlminb","constrOptim"))){
      param <- hyperbChangePars(1, 2,
                  param = c(param[1], exp(param[2]), param[3], exp(param[4])))
    } else {
      param <- hyperbChangePars(1, 2, param = param)
    }

    names(param) <- c("mu", "delta", "alpha", "beta")

    maxLik <- -(as.numeric(opOut[[ind[2]]]))        # maximum likelihood
    conv <- as.numeric(opOut[[ind[4]]])             # convergence
    iter <- as.numeric(opOut[[ind[3]]])[1]          # iterations

  }

  ## Change paramStart back to the primary parameter set version normal scale
  if (!(method %in% c("L-BFGS-B","nlminb","constrOptim"))){
    paramStart <- hyperbChangePars(1, 2,
                    param = c(paramStart[1], exp(paramStart[2]),
                              paramStart[3], exp(paramStart[4])))
  } else {
    paramStart <- hyperbChangePars(1, 2, param = paramStart)
  }

  fitResults <- list(param = param, maxLik = maxLik, criterion = criterion,
                     method = method, conv = conv, iter = iter,
                     obs = x, obsName = xName, paramStart = paramStart,
                     svName = svName, startValues = startValues,
                     breaks = breaks, midpoints = midpoints,
                     empDens = empDens)

  class(fitResults) <- c("hyperbFit", "distFit")

  if (printOut)
    print(fitResults, ...)

  if (plots)
    plot.hyperbFit(fitResults, ...)

  return(fitResults)
} ## End of hyperbFit()


### Function to print object of class hyperbFit
### CYD 01/04/10
### DJS 11/08/06
print.hyperbFit <-
  function(x, digits = max(3, getOption("digits") - 3), ...) {

  if (! "hyperbFit" %in% class(x))
    stop("Object must belong to class hyperbFit")

  cat("\nData:     ", x$obsName, "\n")
  cat("Parameter estimates:\n")
  print.default(format(x$param, digits = digits),
                print.gap = 2, quote = FALSE)
  cat("Likelihood:        ", x$maxLik, "\n")
  cat("criterion :        ", x$criterion , "\n")
  cat("Method:            ", x$method, "\n")
  cat("Convergence code:  ", x$conv, "\n")
  cat("Iterations:        ", x$iter, "\n")
  invisible(x)
} ## End of print.hyperbFit

### Function to plot results of fitting a hyperbolic distribution
plot.hyperbFit <- function(x, which = 1:4,
                           plotTitles = paste(c("Histogram of ",
                                                "Log-Histogram of ",
                                                "Q-Q Plot of ",
                                                "P-P Plot of "),
                                              x$obsName, sep = ""),
                           ask = prod(par("mfcol")) < length(which) &
                                 dev.interactive(), ...) {

  if (! "hyperbFit" %in% class(x))
    stop("Object must belong to class hyperbFit")

  if (ask) {
    op <- par(ask = TRUE)
    on.exit(par(op))
  }

  par(mar = c(6, 4, 4, 2) + 0.1)
  show <- rep(FALSE, 4)
  show[which] <- TRUE
  param <- x$param
  breaks <- x$breaks
  empDens <- x$empDens
  mipoints <- x$midpoints
  obs <- x$obs
  obsName <- x$obsName

  hypDens <- function(x)
    dhyperb(x, param = param)

  logHypDens <- function(x)
    log(dhyperb(x, param = param))

  ymax <- 1.06 * max(hypDens(seq(min(breaks), max(breaks), 0.1)),
                     empDens, na.rm = TRUE)
  if (show[1]) {
    hist.default(obs, breaks, right = FALSE, freq = FALSE, ylim = c(0, ymax),
                 main = plotTitles[1], ...)
    curve(hypDens, min(breaks) - 1, max(breaks) + 1, add = TRUE, ylab = NULL)
    title(sub = paste("param = (",
          round(param[1], 3), ", ", round(param[2], 3), ", ",
          round(param[3], 3), ", ", round(param[4], 3), ")", sep = ""))
  }

  if (show[2]) {
    logHist(obs, breaks, include.lowest = TRUE, right = FALSE,
            main = plotTitles[2], ...)
    curve(logHypDens, min(breaks) - 1, max(breaks) + 1, add = TRUE,
          ylab = NULL, xlab = NULL)
    title(sub = paste("param = (",
          round(param[1], 3), ", ", round(param[2], 3), ", ",
          round(param[3], 3), ", ", round(param[4], 3), ")", sep = ""))
  }

  if (show[3])
    qqhyperb(obs, param = param, main = plotTitles[3], ...)

  if (show[4])
    pphyperb(obs, param = param, main = plotTitles[4], ...)

  invisible()
}

coef.hyperbFit <- function(object, ...) {
  object$param
}

vcov.hyperbFit <- function(object, ...) {
  obs <- object$obs
  param <- object$param
  hessian <- hyperbHessian(obs, param, hessianMethod= "exact",
                           whichParam = 2)
  varcov <- solve(hessian)
  varcov
}

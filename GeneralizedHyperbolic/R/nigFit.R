### Function to fit normal inverse Gaussian distribution to data
### CYD 07/06/10
### DJS 11/09/06
nigFit <- function(x, freq = NULL, paramStart = NULL,
                   startMethod = c("Nelder-Mead","BFGS"),
                   startValues = c("FN","Cauchy","MoM","US"),
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
  startInfo <- nigFitStart(x, startValues = startValues,
                           paramStart = paramStart,
                           startMethodMoM = startMethod, ...)
  paramStart <- startInfo$paramStart
  ## change paramStart in the log scale of param set number 1 (mu,delta,pi,zeta)
  paramStart <- hyperbChangePars(2, 1, param = paramStart)
  if (!(method %in% c("L-BFGS-B","nlminb","constrOptim"))){
    paramStart <- c(paramStart[1], log(paramStart[2]),
                    paramStart[3], log(paramStart[4]))
  }
  svName <- startInfo$svName
  breaks <- startInfo$breaks
  empDens <- startInfo$empDens
  midpoints <- startInfo$midpoints


  if (criterion == "MLE") {
    if (!(method %in% c("L-BFGS-B","nlminb","constrOptim"))){
      llfunc <- function(param) {
        mu <- param[1]
        delta <- exp(param[2])
        hyperbPi <- param[3]
        zeta <- exp(param[4])
        KNu <- besselK(zeta, nu = -1/2)
        KNu2 <- besselK((zeta*sqrt(1 + hyperbPi^2)/delta)*
                sqrt(delta^2 + (x - mu)^2), nu = -1)
        nigDens <- sqrt(1 + hyperbPi^2)*zeta^(1/2)*
                   (sqrt(2*pi)*KNu)^(-1)*(delta^2 + (x - mu)^2)^(-1/2)*
                   KNu2*exp((zeta*hyperbPi/delta)*(x - mu))
        return(-sum(log(nigDens)))
      }
    } else {
      llfunc <- function(param) {
        mu <- param[1]
        delta <- param[2]
        hyperbPi <- param[3]
        zeta <- param[4]
        KNu <- besselK(zeta, nu = -1/2)
        KNu2 <- besselK((zeta*sqrt(1 + hyperbPi^2)/delta)*
                sqrt(delta^2 + (x - mu)^2), nu = -1)
        nigDens <- sqrt(1 + hyperbPi^2)*zeta^(1/2)*
                   (sqrt(2*pi)*KNu)^(-1)*(delta^2 + (x - mu)^2)^(-1/2)*
                   KNu2*exp((zeta*hyperbPi/delta)*(x - mu))
        return(-sum(log(nigDens)))
      }
    }

    output <- numeric(7)
    ind <- 1:4

    if (method == "BFGS") {
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
      cat("Starting loglikelihood = ", llfunc(paramStart, x), " \n")
      opOut <- optim(par = paramStart, llfunc, NULL, x = x,
                     method = "L-BFGS-B",
                     lower = c(-Inf,10^(-2),-Inf,10^(-2)),
                     control = controlLBFGSB, ...)
    }

    if (method == "nlminb") {
      ind <- c(1, 2, 3)
      cat("paramStart =", paramStart[1],paramStart[2],paramStart[3],
          paramStart[4],"\n")
      cat("Starting loglikelihood = ", llfunc(paramStart, x), " \n")
      opOut <- nlminb(start = paramStart, llfunc, NULL, x = x,
                     lower = c(-Inf,0,-Inf,0),
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
    maxLik <- -(as.numeric(opOut[[ind[2]]]))        # maximum likelihood

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

  class(fitResults) <- c("nigFit", "distFit")

  if (printOut)
    print(fitResults, ...)

  if (plots)
    plot.nigFit(fitResults, ...)

  return(fitResults)
} ## End of nigFit()


### Function to print object of class nigFit
### CYD 01/04/10
### DJS 11/08/06
print.nigFit <-
  function(x, digits = max(3, getOption("digits") - 3), ...) {

  if (! "nigFit" %in% class(x))
    stop("Object must belong to class nigFit")

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
} ## End of print.nigFit

### Function to plot results of fitting a nig distribution
plot.nigFit <- function(x, which = 1:4,
                           plotTitles = paste(c("Histogram of ",
                                                "Log-Histogram of ",
                                                "Q-Q Plot of ",
                                                "P-P Plot of "),
                                              x$obsName, sep = ""),
                           ask = prod(par("mfcol")) < length(which) &
                                 dev.interactive(), ...) {

  if (! "nigFit" %in% class(x))
    stop("Object must belong to class nigFit")

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

  nigDens <- function(x)
    dnig(x, param = param)

  logNigDens <- function(x)
    log(dnig(x, param = param))

  ymax <- 1.06 * max(nigDens(seq(min(breaks), max(breaks), 0.1)),
                     empDens, na.rm = TRUE)
  if (show[1]) {
    hist.default(obs, breaks, right = FALSE, freq = FALSE, ylim = c(0, ymax),
                 main = plotTitles[1], ...)
    curve(nigDens, min(breaks) - 1, max(breaks) + 1, add = TRUE, ylab = NULL)
    title(sub = paste("param = (",
          round(param[1], 3), ", ", round(param[2], 3), ", ",
          round(param[3], 3), ", ", round(param[4], 3), ")", sep = ""))
  }

  if (show[2]) {
    logHist(obs, breaks, include.lowest = TRUE, right = FALSE,
            main = plotTitles[2], ...)
    curve(logNigDens, min(breaks) - 1, max(breaks) + 1, add = TRUE,
          ylab = NULL, xlab = NULL)
    title(sub = paste("param = (",
          round(param[1], 3), ", ", round(param[2], 3), ", ",
          round(param[3], 3), ", ", round(param[4], 3), ")", sep = ""))
  }

  if (show[3])
    qqnig(obs, param = param, main = plotTitles[3], ...)

  if (show[4])
    ppnig(obs, param = param, main = plotTitles[4], ...)

  invisible()
}

coef.nigFit <- function(object, ...) {
  object$param
}

vcov.nigFit <- function(object, ...) {
  obs <- object$obs
  param <- object$param
  hessian <- nigHessian(obs, param, hessianMethod= "exact",
                           whichParam = 2)
  varcov <- solve(hessian)
  varcov
}

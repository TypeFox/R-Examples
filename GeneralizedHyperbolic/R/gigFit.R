### Function to fit the Generalized Inverse Gaussian distribution to data
gigFit <- function(x, freq = NULL, paramStart = NULL,
                   startMethod = c("Nelder-Mead","BFGS"),
                   startValues = c("LM","GammaIG","MoM","Symb","US"),
                   method = c("Nelder-Mead","BFGS","nlm"),
                   stand = TRUE, plots = FALSE, printOut = FALSE,
                   controlBFGS = list(maxit = 200),
                   controlNM = list(maxit = 1000),
                   maxitNLM = 1500, ...) {

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


  ## Standardize data?
  sx <- 1
  if (stand) {
    sx <- sd(x)
  }
  x <- x/sx

  startInfo <- gigFitStart(x, startValues = startValues,
                           paramStart = paramStart,
                           startMethodMoM = startMethod, ...)
  paramStart <- startInfo$paramStart
  svName <- startInfo$svName
  breaks <- startInfo$breaks
  empDens <- startInfo$empDens
  midpoints <- startInfo$midpoints

 ## Change paramStart into the log scale

  paramStart <- c(log(paramStart[1]), log(paramStart[2]), paramStart[3])

  ## Create the Log Likelihood function

  llfunc <- function(param) {
    loggigDens <- param[3]/2*log(exp(param[2])/exp(param[1])) -
      log(2*besselK(sqrt(exp(param[1])*exp(param[2])), nu = param[3])) +
        (param[3] - 1)*log(x) - 1/2*(exp(param[1])*x^-1 + exp(param[2])*x)
    as.numeric(loggigDens)
    return(-sum(loggigDens))
  }

  ind <- 1:4

  if (method == "Nelder-Mead") {
    opOut <- optim(paramStart, llfunc, NULL, method = "Nelder-Mead",
                   control = controlNM, ...)
  }

  if (method == "nlm") {
    ind <- c(2, 1)
    opOut <- nlm(llfunc, paramStart, iterlim = maxitNLM, ...)
  }

  param <- as.numeric(opOut[[ind[1]]])[1:3]
  if (stand) {
      param <- c(sx*exp(param[1]), exp(param[2])/sx, param[3])
  } else {
      param <- c(exp(param[1]), exp(param[2]), param[3])
  }

  names(param) <- c("chi", "psi", "lambda")
  maxLik <- -(as.numeric(opOut[[ind[2]]]))
  conv <- as.numeric(opOut[[ind[4]]])
  iter <- as.numeric(opOut[[ind[3]]])[1]

  fitResults <- list(param = param, maxLik = maxLik,
                     method = method, conv = conv,
                     iter = iter, obs = x*sx,
                     obsName = xName, paramStart = paramStart,
                     svName = svName, startValues = startValues,
                     breaks = breaks, midpoints = midpoints,
                     empDens = empDens)
  class(fitResults) <- c("gigFit", "distFit")
  if (printOut)
    print(fitResults, ...)
  if (plots)
    plot.gigFit(fitResults, ...)
  return(fitResults)
}

### print method for gigFit
print.gigFit <- function(x, digits = max(3, getOption("digits") - 3), ...)
{

  if (! "gigFit" %in% class(x))
    stop("Object must belong to class gigFit")

  cat("\nData:     ", x$obsName, "\n")
  cat("Parameter estimates:\n")
  print.default(format(x$param, digits = digits),
                print.gap = 2, quote = FALSE)
  cat("Likelihood:        ", x$maxLik, "\n")
  cat("Method:            ", x$method, "\n")
  cat("Convergence code:  ", x$conv, "\n")
  cat("Iterations:        ", x$iter, "\n")
  invisible(x)
}


### plot method for gigFit
plot.gigFit <- function(x, which = 1:4,
                        plotTitles = paste(c("Histogram of ",
                        "Log-Histogram of ",
                        "Q-Q Plot of ",
                        "P-P Plot of "),
                        x$obsName, sep = ""),
                        ask = prod(par("mfcol")) < length(which) &
                        dev.interactive(), ...) {

  if (! "gigFit" %in% class(x))
    stop("Object must belong to class gigFit")

  if (ask) {
    op <- par(ask = TRUE)
    on.exit(par(op))
  }

  par(mar = c(6, 4, 4, 2) + 0.1)
  show <- rep(FALSE, 4)
  show[which] <- TRUE
  param <- x$param
  obs <- x$obs

  if (show[1]) {
    hist.default(obs, right = FALSE, freq = FALSE,
                 main = plotTitles[1], ...)
    curve(dgig(x, param = param), add = TRUE, ylab = NULL)
    title(sub = paste("param = (", round(param[1], 3), ", ",
          round(param[2], 3), ", ", round(param[3], 3), ")", sep = ""))
  }

  if (show[2]) {
    logHist(obs, include.lowest = TRUE, right = FALSE,
            main = plotTitles[2], ...)
    curve(log(dgig(x, param = param)), add = TRUE, ylab = NULL, xlab = NULL)
    title(sub = paste("param = (", round(param[1], 3), ", ",
          round(param[2], 3), ", ", round(param[3], 3), ")", sep = ""))
  }

  if (show[3])
    qqgig(obs, param = param, main = plotTitles[3], ...)

  if (show[4])
    ppgig(obs, param = param, main = plotTitles[4], ...)

  invisible()
}

### coef method for gigFit
coef.gigFit <- function(object, ...) {
  object$param
}

### vcov method for gigFit
vcov.gigFit <- function(object, ...) {
  obs <- object$obs
  param <- object$param
  hessian <- gigHessian(obs, param, hessianMethod= "tsHessian")
  varcov <- solve(hessian)
  varcov
}


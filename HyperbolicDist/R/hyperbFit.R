### Function to fit hyperbolic distribution to data
###
### DJS 11/09/06
hyperbFit <- function(x, freq = NULL, breaks = NULL, ThetaStart = NULL,
                      startMethod = "Nelder-Mead", startValues = "BN",
                      method = "Nelder-Mead", hessian = FALSE,
                      plots = FALSE, printOut = FALSE,
                      controlBFGS = list(maxit=200),
                      controlNM = list(maxit=1000), maxitNLM = 1500, ...)
{
  xName <- paste(deparse(substitute(x), 500), collapse = "\n")
  if(!is.null(freq)){
      if(length(freq)!=length(x)){
      stop("vectors x and freq are not of the same length")
    }
    x <- rep(x, freq)
  }
  x <- as.numeric(na.omit(x))
  startInfo <- hyperbFitStart(x, breaks = breaks,
                              startValues = startValues,
                              ThetaStart = ThetaStart,
                              startMethodSL = startMethod,
                              startMethodMoM = startMethod,...)
  ThetaStart <- startInfo$ThetaStart
  svName <- startInfo$svName
  breaks <- startInfo$breaks
  empDens <- startInfo$empDens
  midpoints <- startInfo$midpoints

  llfunc <- function(Theta){
    KNu <- besselK(exp(Theta[2]), nu = 1)
    -sum(log(dhyperb(x, Theta, KNu = KNu, logPars = TRUE)))
  }
  output <- numeric(7)
  ind <- 1:4
  if(method == "BFGS"){
    opOut <- optim(ThetaStart, llfunc, NULL, method = "BFGS",
                  hessian = hessian, control = controlBFGS,...)
  }

  if(method == "Nelder-Mead"){
    opOut <- optim(ThetaStart, llfunc, NULL, method = "Nelder-Mead",
                  hessian = hessian, control = controlNM,...)
  }

  if(method == "nlm"){
    ind <- c(2,1,5,4)
    opOut <- nlm(llfunc, ThetaStart, hessian = hessian,
                iterlim = maxitNLM,...)

  }

  Theta <- as.numeric(opOut[[ind[1]]])[1:4]       # parameter values
  Theta[2] <- exp(Theta[2])                       # don't use logs
  Theta[3] <- exp(Theta[3])                       # don't use logs
  names(Theta) <- c("pi","zeta","delta","mu")
  maxLik <- -(as.numeric(opOut[[ind[2]]]))        # maximum likelihood
  conv <- as.numeric(opOut[[ind[4]]])             # convergence
  iter <- as.numeric(opOut[[ind[3]]])[1]          # iterations
  ThetaStart <- c(ThetaStart[1], exp(ThetaStart[2]),
                  exp(ThetaStart[3]), ThetaStart[4])

  KNu <- besselK(Theta[2], nu = 1)

  fitResults <- list(Theta = Theta, maxLik = maxLik,
                     hessian = if (hessian) opOut$hessian else NULL,
                     method = method, conv = conv, iter = iter,
                     obs = x, obsName = xName, ThetaStart = ThetaStart,
                     svName = svName, startValues = startValues,
                     KNu = KNu, breaks = breaks,
                     midpoints = midpoints, empDens = empDens)

  class(fitResults) <- "hyperbFit"

  if(printOut == TRUE){
    print(fitResults, ...)
  }
  if(plots == TRUE){
    plot.hyperbFit(fitResults, ...)
  }
  fitResults
} ## End of hyperbFit()


### Function to print object of class hyperbFit
### DJS 11/08/06
print.hyperbFit <- function(x,
                            digits = max(3, getOption("digits") - 3), ...)
{
  if (!class(x) == "hyperbFit"){
    stop("Object must belong to class hyperbFit")
  }
  cat("\nData:     ", x$obsName, "\n")
  cat("Parameter estimates:\n")
  print.default(format(x$Theta, digits = digits),
            print.gap = 2, quote = FALSE)
  cat("Likelihood:        ", x$maxLik, "\n")
  cat("Method:            ", x$method, "\n")
  cat("Convergence code:  ", x$conv, "\n")
  cat("Iterations:        ", x$iter, "\n")
  invisible(x)
} ## End of print.hyperbFit

### Function to plot results of fitting a hyperbolic distribution
plot.hyperbFit <-
function(x, which = 1:4,
         plotTitles=paste(c("Histogram of ","Log-Histogram of ",
                            "Q-Q Plot of ","P-P Plot of "),
                            x$obsName, sep=""),
         ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...){
  if (!class(x)=="hyperbFit")
    stop("Object must belong to class hyperbFit")
  if (ask) {
    op <- par(ask = TRUE)
    on.exit(par(op))
  }
  par(mar=c(6,4,4,2)+0.1)
  show <- rep(FALSE, 4)
  show[which] <- TRUE
  Theta <- x$Theta
  KNu <- x$KNu
  breaks <- x$breaks
  empDens <- x$empDens
  mipoints <- x$midpoints
  obs <- x$obs
  obsName <- x$obsName

  hypDens <- function(x) dhyperb(x, Theta, KNu, logPars = FALSE)

  logHypDens <- function(x) log(dhyperb(x, Theta, KNu = KNu, logPars = FALSE))

  ymax <- 1.06 * max(hypDens(seq(min(breaks), max(breaks), 0.1)),
                     empDens, na.rm = TRUE)
  if (show[1]){
    hist.default(obs, breaks, right = FALSE, freq = FALSE, ylim = c(0, ymax),
                 main = plotTitles[1], ...)
    curve(hypDens, min(breaks) - 1, max(breaks) + 1, add = TRUE, ylab = NULL)
    title(sub = paste("Theta = (",
          round(Theta[1], 3), ",", round(Theta[2], 3), ",",
          round(Theta[3], 3), ",", round(Theta[4], 3), ")", sep = ""))
  }
  if (show[2]){
    logHist(obs, breaks, include.lowest = TRUE, right = FALSE,
            main = plotTitles[2], ...)
    curve(logHypDens, min(breaks) - 1, max(breaks) + 1, add = TRUE,
          ylab = NULL, xlab = NULL)
    title(sub = paste("Theta = (",
          round(Theta[1], 3), ",", round(Theta[2], 3), ",",
          round(Theta[3], 3), ",", round(Theta[4], 3), ")", sep = ""))
  }

  if (show[3]){
    qqhyperb(obs, Theta, main = plotTitles[3], ...)
  }
  if (show[4]){
    pphyperb(obs, Theta, main = plotTitles[4], ...)
  }
  invisible()
}

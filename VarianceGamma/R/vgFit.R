vgFit <- function (x, freq = NULL, breaks = NULL, paramStart = NULL,
                   startMethod = "Nelder-Mead", startValues = "SL",
                   method = "Nelder-Mead", hessian = FALSE,
                   plots = FALSE, printOut = FALSE,
                   controlBFGS = list(maxit = 200),
                   controlNM = list(maxit = 1000),
                   maxitNLM = 1500, ...) {
    xName <- paste(deparse(substitute(x), 500), collapse = "\n")
    if (!is.null(freq)) {
        if (length(freq) != length(x)) {
            stop("vectors x and freq are not of the same length")
        }
        x <- rep(x, freq)
    }
    x <- as.numeric(na.omit(x))
    startInfo <- vgFitStart(x = x, breaks = breaks, startValues = startValues,
                            paramStart = paramStart,
                            startMethodSL = startMethod,
                            startMethodMoM = startMethod, ...)
    paramStart <- startInfo$paramStart
    svName <- startInfo$svName
    breaks <- startInfo$breaks
    empDens <- startInfo$empDens
    midpoints <- startInfo$midpoints
    llfunc <- function(logParam) {
      vgC <- logParam[1]
      sigma <- exp(logParam[2])
      theta <- logParam[3]
      nu <- exp(logParam[4])
      ## added by Christine, email 23/09/2009
      adjust <- 0.001
      if (sigma < 0 | (abs(sigma - 0) < adjust)) {
        sigma <- adjust
      }
      if (nu < 0 | (abs(nu - 0) < adjust)) {
        nu <- adjust
      }
      ## end of addition
      return(-sum(log(dvg(x = x, param = c(vgC, sigma, theta, nu)))))
    }
    output <- numeric(7)
    ind <- 1:4
    if (method == "BFGS") {
        opOut <- optim(paramStart, llfunc, NULL, method = "BFGS",
                       hessian = hessian, control = controlBFGS, ...)
    }
    if (method == "Nelder-Mead") {
        opOut <- optim(paramStart, llfunc, NULL, method = "Nelder-Mead",
                       hessian = hessian, control = controlNM, ...)
    }
    if (method == "nlm") {
        ind <- c(2, 1, 5, 4)
        opOut <- nlm(llfunc, paramStart, hessian = hessian,
                     iterlim = maxitNLM, ...)
    }
    param <- as.numeric(opOut[[ind[1]]])[1:4]
    param[2] <- exp(param[2])
    param[4] <- exp(param[4])
    names(param) <- c("vgC", "sigma", "theta", "nu")
    maxLik <- -(as.numeric(opOut[[ind[2]]]))
    conv <- as.numeric(opOut[[ind[4]]])
    iter <- as.numeric(opOut[[ind[3]]])[1]
    paramStart <- c(paramStart[1], exp(paramStart[2]),
                    paramStart[3], exp(paramStart[4]))
    fitResults <- list(param = param, maxLik = maxLik,
                       hessian = if (hessian) opOut$hessian else NULL,
                       method = method, conv = conv, iter = iter, obs = x,
                       obsName = xName, paramStart = paramStart,
                       svName = svName, startValues = startValues,
                       breaks = breaks, midpoints = midpoints,
                       empDens = empDens)
    class(fitResults) <- "vgFit"
    if (printOut == TRUE) {
        print(fitResults, ...)
    }
    if (plots == TRUE) {
        plot.vgFit(fitResults, ...)
    }
    return(fitResults)
}

print.vgFit <- function (x, digits = max(3, getOption("digits") - 3), ...)
{
  if (!class(x) == "vgFit") {
    stop("Object must belong to class vgFit")
  }
  cat("\nData:     ", x$obsName, "\n")
  cat("Parameter estimates:\n")
  print.default(format(x$param, digits = digits), print.gap = 2,
                quote = FALSE)
  cat("Likelihood:        ", x$maxLik, "\n")
  cat("Method:            ", x$method, "\n")
  cat("Convergence code:  ", x$conv, "\n")
  cat("Iterations:        ", x$iter, "\n")
  invisible(x)
}

plot.vgFit <- function (x, which = 1:4,
                        plotTitles = paste(c("Histogram of ",
                        "Log-Histogram of ", "Q-Q Plot of ", "P-P Plot of "),
                        x$obsName, sep = ""),
                        ask = prod(par("mfcol")) < length(which) &&
                        dev.interactive(), ...)
{
    if (!class(x) == "vgFit")
        stop("Object must belong to class vgFit")
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
    vgDens <- function(obs) dvg(obs, param = param,
                                tolerance = .Machine$double.eps ^ 0.5)
    logvgDens <- function(obs) log(dvg(obs, param = param,
                                       tolerance = .Machine$double.eps ^ 0.5))
    ymax <- 1.06 * max(vgDens(seq(min(breaks), max(breaks), 0.1)),
                       empDens, na.rm = TRUE)
    if (show[1]) {
        hist.default(obs, breaks, right = FALSE, freq = FALSE,
                     ylim = c(0, ymax), main = plotTitles[1], ...)
        curve(vgDens, min(breaks) - 1, max(breaks) + 1, add = TRUE,
              ylab = NULL)
        title(sub = paste("param = (", round(param[1], 3), ",",
              round(param[2], 3), ",", round(param[3], 3), ",",
              round(param[4], 3), ")", sep = ""))
    }
    if (show[2]) {
        logHist(obs, breaks, include.lowest = TRUE, right = FALSE,
                main = plotTitles[2], ...)
        curve(logvgDens, min(breaks) - 1, max(breaks) + 1, add = TRUE,
              ylab = NULL, xlab = NULL)
        title(sub = paste("param = (", round(param[1], 3), ",",
              round(param[2], 3), ",", round(param[3], 3), ",",
              round(param[4], 3), ")", sep = ""))
    }
    if (show[3]) {
      qqvg(obs, param = param, main = plotTitles[3], ...)
    }
    if (show[4]) {
      ppvg(obs, param = param, main = plotTitles[4], ...)
    }
    invisible()
}



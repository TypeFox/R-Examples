qqvg <- function (y, vgC = NULL, sigma = NULL, theta = NULL, nu = NULL,
                  param = c(vgC,sigma,theta,nu),
                  main = "Variance Gamma Q-Q Plot",
                  xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
                  plot.it = TRUE, line = TRUE, ...) {

  if (has.na <- any(ina <- is.na(y))) {
    yN <- y
    y <- y[!ina]
  }
  if (0 == (n <- length(y))){
    stop("y is empty or has only NAs")
  }

  if (length(param) > 0) {
    #check parameters
    parResult <- vgCheckPars(param = param)
    case <- parResult$case
    errMessage <- parResult$errMessage
    if (case == "error"){
      stop(errMessage)
    }
    vgC <- param[1]
    sigma <- param[2]
    theta <- param[3]
    nu <- param[4]
  } else {
    fitResults <- vgFit(y, freq = NULL, breaks = NULL, paramStart = NULL,
                        startMethod = "Nelder-Mead", startValues = "SL",
                        method = "Nelder-Mead", hessian = FALSE,
                        plots = FALSE, printOut = FALSE,
                        controlBFGS = list(maxit = 200),
                        controlNM = list(maxit = 1000),
                        maxitNLM = 1500, ...)
    param <- fitResults$param
    names(param) = NULL
    vgC <- param[1]
    sigma <- param[2]
    theta <- param[3]
    nu <- param[4]
  }

  x <- qvg(ppoints(n), param = param)[order(order(y))]

  if (has.na) {
    y <- x
    x <- yN
    x[!ina] <- y
    y <- yN
  }
  if (plot.it) {
    qqplot(x, y, main = main, xlab = xlab, ylab = ylab, ...)
    title(sub = paste("vgC = ", round(vgC, 3), ", sigma = ",
          round(sigma, 3), ", theta = ", round(theta, 3),
          ", nu = ", round(nu, 3), sep = ""))
  }
  if (line) abline(0, 1)
  invisible(list(x = x, y = y))
}

ppvg <- function (y, vgC = NULL, sigma = NULL, theta = NULL, nu = NULL,
                  param = c(vgC,sigma,theta,nu),
                  main = "Variance Gamma P-P Plot",
                  xlab = "Uniform Quantiles",
                  ylab = "Probability-integral-transformed Data",
                  plot.it = TRUE, line = TRUE, ...) {
  if (has.na <- any(ina <- is.na(y))) {
    yN <- y
    y <- y[!ina]
  }
  if (0 == (n <- length(y))) {
    stop("data is empty")
  }

  if (length(param) > 0) {
    #check parameters
    parResult <- vgCheckPars(param = param)
    case <- parResult$case
    errMessage <- parResult$errMessage
    if (case == "error"){
      stop(errMessage)
    }
    vgC <- param[1]
    sigma <- param[2]
    theta <- param[3]
    nu <- param[4]
  } else {
    fitResults <- vgFit(y, freq = NULL, breaks = NULL, paramStart = NULL,
                        startMethod = "Nelder-Mead", startValues = "SL",
                        method = "Nelder-Mead", hessian = FALSE,
                        plots = FALSE, printOut = FALSE,
                        controlBFGS = list(maxit = 200),
                        controlNM = list(maxit = 1000),
                        maxitNLM = 1500, ...)
    param <- fitResults$param
    names(param) = NULL
    vgC <- param[1]
    sigma <- param[2]
    theta <- param[3]
    nu <- param[4]
  }

  yvals <- pvg(y, param = param)
  xvals <- ppoints(n, a = 1/2)[order(order(y))]
  if (has.na) {
    y <- yvals
    x <- xvals
    yvals <- yN
    yvals[!ina] <- y
    xvals <- yN
    xvals[!ina] <- x
  }
  if (plot.it) {
    plot(xvals, yvals, main = main, xlab = xlab, ylab = ylab,
         ylim = c(0, 1), xlim = c(0, 1), ...)
    title(sub = paste("vgC=", round(vgC, 3), ", sigma=",
          round(sigma, 3), ", theta=", round(theta, 3),
          ", nu=", round(nu, 3), sep = ""))
  }
  if (line) abline(0, 1)
  invisible(list(x = xvals, y = yvals))
}




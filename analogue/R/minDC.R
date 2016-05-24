###########################################################################
##                                                                       ##
## minDC & methods - function to extract minimum DC values from objects  ##
##                   and assocaited methods                              ##
##                                                                       ##
## Created       : 12-Feb-2007                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.0-1                                                 ##
## Last modified : 12-Feb-2007                                           ##
##                                                                       ##
###########################################################################
minDC <- function(x, ...) UseMethod("minDC")

minDC.default <- function(x, ...) {
  if(is.null(x$minDC))
    stop("No 'minDC' component of 'x' found.")
  res <- list(minDC = x$minDC)
  class(res) <- "minDC"
  return(res)
}

minDC.predict.mat <- function(x, ...) {
  if(is.null(x$minDC))
    stop("'minDC' is not currently available for all\n'predict.mat' options.")
  res <- list(minDC = x$minDC, method = x$method, quantiles = x$quantiles)
  class(res) <- "minDC"
  return(res)
}

minDC.analog <- function(x, probs = c(0.01, 0.02, 0.05, 0.1),
                         ...) {
  if(class(x) != "analog")
    stop("'x' is not of class \"analog\".")
  minDC <- apply(x$analogs, 2, min)
  quant <- NULL
  if(!is.null(x$train))
    quant <- quantile(as.dist(x$train), probs = probs)
  res <- list(minDC = minDC, method = x$method, quantiles = quant)
  class(res) <- "minDC"
  return(res)
}

print.minDC <- function(x, digits = min(3, getOption("digits") - 4), ...) {
  cat("\n")
  writeLines(strwrap("Minimum dissimilarity per sample", prefix = "\t"))
  cat("\n")
  cat(paste("Dissimilarity:", x$method, "\n\n"))
  print(unclass(x$minDC), digits = digits)
}

plot.minDC <- function(x, depths, use.labels = FALSE,
                       quantiles = TRUE,
                       rev.x = TRUE,
                       type = "l",
                       xlim, ylim,
                       xlab = "", ylab = "Dissimilarity",
                       main = "",
                       sub = NULL,
                       col.quantile = "red",
                       lty.quantile = "dotted",
                       ...) {
  if(quantiles) {
    mars <- par("mar")
    opar <- par(mar = c(mars[1:3], 3.1))
    on.exit(par(opar))
  }
  if(missing(depths))
    {
      if(use.labels) {
        depths <- as.numeric(names(x$minDC))
      } else {
        stop("If \"use.labels = FALSE\", then \"depths\" must be provided.")
      }
    }
  if(missing(xlim))
    xlim <- range(depths)
  if(rev.x)
    xlim <- rev(xlim)
  if(missing(ylim))
    ylim <- range(x$minDC)
  if(is.null(sub))
    sub <- paste("Dissimilarity:", x$method)
  plot(depths, x$minDC, ylim = ylim, xlim = xlim, type = "n",
       ylab = ylab, xlab = xlab, main = main, sub = sub, ...)
  if(quantiles & !is.null(x$quantile)) {
    abline(h = x$quantile, lty = lty.quantile, col = col.quantile)
    axis(4, at = x$quantile, labels = names(x$quantile), las = 2)
  }
  lines(depths, x$minDC, type = type)
  invisible()
}

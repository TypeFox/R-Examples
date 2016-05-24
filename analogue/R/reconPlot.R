###########################################################################
##                                                                       ##
## reconPlot - function to draw reconstructed environmental variables    ##
##             from transfer function models                             ##
##                                                                       ##
## Created       : 05-Nov-2006                                           ##
## Author        : Gavin Simpson                                         ##
## Version       : 0.1-0                                                 ##
## Last modified : 03-Mar-2007                                           ##
##                                                                       ##
###########################################################################
##
## Arguments
##
## S3 method
reconPlot <- function(x, ...) UseMethod("reconPlot")

reconPlot.default <- function(x, depths, errors,
                              display.error = c("none", "bars", "lines"),
                              rev.x = TRUE,
                              col.error = "grey", lty.error = "dashed",
                              type = "l",
                              xlim, ylim,
                              xlab = "", ylab = "", main = "",
                              ...) {
  ##stop("No default method for \"reconPlot\"")
  if(missing(display.error))
    display.error <- "none"
  display.error <- match.arg(display.error)
  show.errors <- FALSE
  if(display.error != "none")
    show.errors <- TRUE
  if(missing(errors) & show.errors)
    stop("'errors' must be supplied if 'display.error != \"none\".")
  if(missing(xlim))
    xlim <- range(depths)
  if(rev.x)
    xlim <- rev(xlim)
  if(show.errors) {
    upper <- x + errors
    lower <- x - errors
  }
  if(missing(ylim)) {
    if(show.errors)
      ylim <- range(x, upper, lower)
    else
      ylim <- range(x)
  }
  plot(depths, x, ylim = ylim, xlim = xlim, type = "n",
       ylab = ylab, xlab = xlab, main = main, ...)
  if(show.errors){
    if(display.error == "bars")
      arrows(depths, upper, depths, lower, length = 0.02, angle = 90,
             code = 3, col = col.error)
    else {
      lines(depths, upper, type = type, lty = lty.error, col = col.error)
      lines(depths, lower, type = type, lty = lty.error, col = col.error)
    }
  }
  lines(depths, x, type = type, ...)
  invisible()
}

reconPlot.predict.mat <- function(x, depths, use.labels = FALSE,
                                  predictions = c("model",
                                    "bootstrap"),
                                  display.error = c("none", "bars", "lines"),
                                  sample.specific = TRUE,
                                  ...) {
    if(missing(display.error))
        display.error <- "none"
    display.error <- match.arg(display.error)
    if(missing(predictions))
        predictions <- "model"
    predictions <- match.arg(predictions)
    if(missing(depths)) {
        if(use.labels)
            depths <- as.numeric(colnames(x$predictions$model$predicted))
        else
            stop("If \"use.labels = FALSE\", then \"depths\" must be provided.")
    }
    if(predictions == "model") {
        n.analogues <- x$predictions$model$k
        preds <- x$predictions$model$predicted[n.analogues, ]
        errors <- x$model$rmsep[n.analogues]
    } else {
        n.analogues <- x$predictions$bootstrap$k
        preds <- x$predictions$bootstrap$predicted[,n.analogues]
        if(sample.specific)
            errors <- x$predictions$sample.errors$rmsep[, n.analogues]
        else
            errors <- x$bootstrap$rmsep[n.analogues]
    }
    reconPlot.default(preds, depths, errors, display.error = display.error, ...)
}

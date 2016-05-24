#' Plotting Confidence/Prediction Bands
#' 
#' Plots fitted model for an object of class \code{"lm"} or \code{"nls"} with 
#' the option of adding a confidence and/or prediction band. 
#'
#' @param object An object that inherits from class \code{"lm"}, \code{"glm"}, 
#'               or \code{"nls"}.
#' @param type The type of prediction required. The default is on the scale of 
#'             the response variable; the alternative \code{"link"} is on the 
#'             scale of the linear predictor. This option is only used when
#'             plotting \code{"glm"} objects.
#' @param interval A character string indicating if a prediction band, 
#'   confidence band, both, or none should be plotted.
#' @param level The desired confidence level.
#' @param data An optional data frame containing the variables in the model. 
#' @param adjust A character string indicating the type of adjustment (if any) 
#' to make to the confidence/prediction bands.
#' @param k An integer to be used in computing the critical value for the 
#' confidence/prediction bands. Only needed when \code{adjust = "Bonferroni"} or
#' when \code{adjust = "Scheffe"} and \code{interval = "prediction"}.
#' @param shade A logical value indicating if the band should be shaded.
#' @param extend.range A logical value indicating if the fitted regression line
#' and bands (if any) should extend to the edges of the plot. Default is 
#' \code{FALSE}.
#' @param col.conf Shade color for confidence band.
#' @param col.pred Shade color for prediction band.
#' @param col.fit The color to use for the fitted line.
#' @param border.conf The color to use for the confidence band border.
#' @param border.pred The color to use for the prediction band border. 
#' @param lty.conf Line type to use for confidence band border.
#' @param lty.pred Line type to use for prediction band border.
#' @param lty.fit Line type to use for the fitted regression line.
#' @param lwd.conf Line width to use for confidence band border.
#' @param lwd.pred Line width to use for prediction band border.
#' @param lwd.fit Line width to use for the fitted regression line.
#' @param n The number of predictor values at which to evaluate the fitted model
#' (larger implies a smoother plot).
#' @param xlab A title for the x axis.
#' @param ylab A title for the y axis.
#' @param xlim The x limits (x1, x2) of the plot.
#' @param ylim The y limits (y1, y2) of the plot. 
#' @param hide A logical value indicating if the fitted model should be 
#' plotted on top of the points (FALSE) or behind them (TRUE). Default is 
#' TRUE.
#' @param ... Additional optional arguments passed on to \code{plot}.
#' @rdname plotFit
#' @importFrom graphics lines plot polygon
#' @importFrom grDevices extendrange grey
#' @importFrom stats family formula getCall predict qnorm
#' @export
#' @note
#' By default, the plotted intervals are pointwise intervals. For simultaneous 
#' intervals use \code{adjust = "Bonferroni"} or \code{adjust = "Scheffe"}. For
#' the Bonferroni adjustment, you must specify a value for \code{k}, the number
#' of intervals for which the coverage is to hold simultaneously. For the 
#' Scheffe adjustment, specifying a value for \code{k} is only required when
#' \code{interval = "prediction"}; if \code{interval = "confidence"}, \code{k} 
#' is set equal to \eqn{p}, the number of regression parameters. For example,
#' if \code{object} is a simple linear regression model, then calling 
#' \code{plotFit} with \code{interval = "confidence"} and 
#' \code{adjust = "Scheffe"} will plot the Working-Hotelling band.
#' 
#' Confidence/prediction bands for nonlinear regression (i.e., objects of class
#' \code{nls}) are based on a linear approximation as described in Bates & Watts 
#' (2007). This funtion was inpired by the \code{\link[nlstools]{plotfit}} function
#' from the \code{nlstools} package.
#' @references
#' Bates, D. M., and Watts, D. G. (2007)
#' \emph{Nonlinear Regression Analysis and its Applications}. Wiley.
#' 
#' F. Baty and M. L. Delignette-Muller (2012), 
#' A Toolbox for Nonlinear Regression in R: The Package nlstools.
#' \emph{Journal of Statistical Software} \bold{(under revision)}.
#' @examples
#' #
#' # A nonlinear regression example
#' #
#' data(Puromycin, package = "datasets")
#' Puromycin2 <- Puromycin[Puromycin$state == "treated", ][, 1:2]
#' Puro.nls <- nls(rate ~ Vm * conc/(K + conc), data = Puromycin2,
#'                 start = c(Vm = 200, K = 0.05))
#' plotFit(Puro.nls, interval = "both", pch = 19, shade = TRUE, 
#'         col.conf = "skyblue4", col.pred = "lightskyblue2")  
plotFit <- function(object, ...) {
  UseMethod("plotFit")
} 


#' @rdname plotFit
#' @method plotFit lm
#' @importFrom graphics plot
#' @importFrom stats formula
#' @export
plotFit.lm <- function(object, 
                       interval = c("none", "both", "confidence", "prediction"), 
                       level = 0.95, data,
                       adjust = c("none", "Bonferroni", "Scheffe"), k, ...,
                       shade = FALSE, extend.range = FALSE, hide = TRUE,
                       col.conf = if (shade) grey(0.7) else "black", 
                       col.pred = if (shade) grey(0.9) else "black",  
                       border.conf = col.conf, border.pred = col.pred, 
                       col.fit = "black", lty.conf = if (shade) 1 else 2, 
                       lty.pred = if (shade) 1 else 3, lty.fit = 1, 
                       lwd.conf = 1, lwd.pred = 1, lwd.fit = 1, n = 500, 
                       xlab, ylab, xlim, ylim) {
  
  # Extract data 
  if (!missing(data)) { 
    .data <- data 
  } else {
    .data <- eval(getCall(object)$data, envir = parent.frame())
  }
  if (is.null(.data)) {  # throw error if no data are found
    stop(paste("No data to plot. If the", class(object), "object does not",
               "contain a data component then a data frame containing the", 
               "variables in the model must be supplied through the 'data'", 
               "argument to plotFit."))
  }
  
  # Extract variable names and values
  xname <- intersect(all.vars(formula(object)[[3]]), colnames(.data)) 
  yname <- all.vars(formula(object)[[2]])
  if (length(xname) != 1) stop("Only one independent variable allowed.")
  if (length(yname) != 1) stop("Only one dependent variable allowed.")
  xvals <- .data[, xname]
  yvals <- with(.data, eval(formula(object)[[2]]))
  
  # Plot limits, labels, etc.
  if (missing(xlim)) xlim <- range(xvals)  # default limits for x-axis
  xgrid <- if (extend.range) {  # the x values at which to evaluate
    list(seq(from = extendrange(xlim)[1], to = extendrange(xlim)[2], 
             length = n))
  } else {
    list(seq(from = xlim[1], to = xlim[2], length = n))
  }
  names(xgrid) <- xname
  if (missing(xlab)) xlab <- xname  # default label for x-axis
  if (missing(ylab)) ylab <- yname  # default label for y-axis

  # Maximum and minimum of fitted values
  interval = match.arg(interval)
  if (interval == "none") {
    fitvals <- predFit(object, newdata = xgrid)
    fit.ymin <- min(fitvals)
    fit.ymax <- max(fitvals)
  }
  
  # Confidence interval for mean response
  adjust <- match.arg(adjust)
  if (interval == "confidence" || interval == "both") {
    conf <- predFit(object, newdata = xgrid, interval = "confidence", 
                     level = level, adjust = adjust, k = k)
    conf.lwr <- conf[, "lwr"]
    conf.upr <- conf[, "upr"]
    conf.ymin <- min(conf.lwr)
    conf.ymax <- max(conf.upr)
  }
  
  # Prediction interval for individual response
  if (interval == "prediction" || interval == "both") {
    pred <- predFit(object, newdata = xgrid, interval = "prediction", 
                     level = level, adjust = adjust, k = k)
    pred.lwr <- pred[, "lwr"]
    pred.upr <- pred[, "upr"]
    pred.ymin <- min(pred.lwr)
    pred.ymax <- max(pred.upr)
  }
  
  # Automatic limits for y-axis
  if (missing(ylim)) {
    if(interval == "prediction" || interval == "both") {
      ylim <- c(min(c(pred.ymin, yvals)), max(c(pred.ymax, yvals)))
    }
    if (interval == "confidence") {
      ylim <- c(min(c(conf.ymin, yvals)), max(c(conf.ymax, yvals)))
    } 
    if (interval == "none") {
      ylim <- c(min(c(fit.ymin, yvals)), max(c(fit.ymax, yvals)))
    }
  }
  
  # Plot data, mean response, etc.
  plot(xvals, yvals, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
       panel.first = if (hide) {  # draw points last
         if (shade) {
           # Draw (hidden) shaded prediction band
           if (interval == "prediction" || interval == "both") {
             polygon(c(xgrid[[1]], rev(xgrid[[1]])), c(pred.lwr, rev(pred.upr)), 
                     col = col.pred, border = border.pred, lty = lty.pred, 
                     lwd = lwd.pred)
           }
           # Draw (hidden) shaded confidence band
           if (interval == "confidence" || interval == "both") {
             polygon(c(xgrid[[1]], rev(xgrid[[1]])), c(conf.lwr, rev(conf.upr)), 
                     col = col.conf, border = border.conf, lty = lty.conf, 
                     lwd = lwd.conf)
           }
         } else {
           # Draw (hidden) unshaded prediction band
           if (interval == "prediction" || interval == "both") {
             lines(xgrid[[1]], pred.lwr, col = col.pred, lty = lty.pred, 
                   lwd = lwd.pred)
             lines(xgrid[[1]], pred.upr, col = col.pred, lty = lty.pred, 
                   lwd = lwd.pred)
           }
           # Draw (hidden) unshaded confidence band
           if (interval == "confidence" || interval == "both") {
             lines(xgrid[[1]], conf.lwr, col = col.conf, lty = lty.conf, 
                   lwd = lwd.conf)
             lines(xgrid[[1]], conf.upr, col = col.conf, lty = lty.conf, 
                   lwd = lwd.conf)
           }
         }
         # Draw (hidden) fitted response curve
         lines(xgrid[[1]], suppressWarnings(predict(object, newdata = xgrid)), 
               lty = lty.fit, lwd = lwd.fit, col = col.fit)  
       } else NULL,
       panel.last = if(!hide) {  # draw points first
         if (shade) {
           # Draw shaded prediction band
           if (interval == "prediction" || interval == "both") {
             polygon(c(xgrid[[1]], rev(xgrid[[1]])), c(pred.lwr, rev(pred.upr)), 
                     col = col.pred, border = border.pred, lty = lty.pred, 
                     lwd = lwd.pred)
           }
           # Draw shaded confidence band
           if (interval == "confidence" || interval == "both") {
             polygon(c(xgrid[[1]], rev(xgrid[[1]])), c(conf.lwr, rev(conf.upr)), 
                     col = col.conf, border = border.conf, lty = lty.conf, 
                     lwd = lwd.conf)
           }
         } else {
           # Draw unshaded prediction band
           if (interval == "prediction" || interval == "both") {
             lines(xgrid[[1]], pred.lwr, col = col.pred, lty = lty.pred, 
                   lwd = lwd.pred)
             lines(xgrid[[1]], pred.upr, col = col.pred, lty = lty.pred, 
                   lwd = lwd.pred)
           }
           # Draw unshaded confidence band
           if (interval == "confidence" || interval == "both") {
             lines(xgrid[[1]], conf.lwr, col = col.conf, lty = lty.conf, 
                   lwd = lwd.conf)
             lines(xgrid[[1]], conf.upr, col = col.conf, lty = lty.conf, 
                   lwd = lwd.conf)
           }
         }
         # Draw fitted response curve
         lines(xgrid[[1]], suppressWarnings(predict(object, newdata = xgrid)), 
               lty = lty.fit, lwd = lwd.fit, col = col.fit)  
       } else NULL, ...)
  
}


#' @rdname plotFit
#' @export
#' @method plotFit nls
plotFit.nls <- function(object, 
                        interval = c("none", "both", "confidence", "prediction"), 
                        level = 0.95, data,
                        adjust = c("none", "Bonferroni", "Scheffe"), k, ..., 
                        shade = FALSE, extend.range = FALSE, hide = TRUE,
                        col.conf = if (shade) grey(0.7) else "black", 
                        col.pred = if (shade) grey(0.9) else "black",  
                        border.conf = col.conf, border.pred = col.pred, 
                        col.fit = "black", lty.conf = if (shade) 1 else 2, 
                        lty.pred = if (shade) 1 else 3, lty.fit = 1, 
                        lwd.conf = 1, lwd.pred = 1, lwd.fit = 1, n = 500, 
                        xlab, ylab, xlim, ylim) {
  
  # Extract data 
  if (!missing(data)) { 
    .data <- data 
  } else {
    .data <- eval(getCall(object)$data, envir = parent.frame())
  }
  if (is.null(.data)) {  # throw error if no data are found
    stop(paste("No data to plot. If the", class(object), "object does not",
               "contain a data component then a data frame containing the", 
               "variables in the model must be supplied through the 'data'", 
               "argument to plotFit."))
  }
  
  # Extract variable names and values
  xname <- intersect(all.vars(formula(object)[[3]]), colnames(.data)) 
  yname <- all.vars(formula(object)[[2]])
  if (length(xname) != 1) stop("Only one independent variable allowed.")
  if (length(yname) != 1) stop("Only one dependent variable allowed.")
  xvals <- .data[, xname]
  yvals <- with(.data, eval(formula(object)[[2]]))
  
  # Plot limits, labels, etc.
  if (missing(xlim)) xlim <- range(xvals)  # default limits for x-axis
  xgrid <- if (extend.range) {  # set up plotting grid
    list(seq(from = extendrange(xlim)[1], to = extendrange(xlim)[2], 
             length = n))
  } else {
    list(seq(from = xlim[1], to = xlim[2], length = n))
  }
  names(xgrid) <- xname
  if (missing(xlab)) xlab <- xname  # default label for x-axis
  if (missing(ylab)) ylab <- yname  # default label for y-axis
  
  # Maximum and minimum of fitted values
  interval = match.arg(interval)
  if (interval == "none") {
    fitvals <- predFit(object, newdata = xgrid)
    fit.ymin <- min(fitvals)
    fit.ymax <- max(fitvals)
  }
  
  # Confidence interval for mean response
  adjust <- match.arg(adjust)
  if (interval == "confidence" || interval == "both") {
    conf <- predFit(object, newdata = xgrid, interval = "confidence", 
                     level = level, adjust = adjust, k = k)
    conf.lwr <- conf[, "lwr"]
    conf.upr <- conf[, "upr"]
    conf.ymin <- min(conf.lwr)
    conf.ymax <- max(conf.upr)
  }
  
  # Prediction interval for individual response
  if (interval == "prediction" || interval == "both") {
    pred <- predFit(object, newdata = xgrid, interval = "prediction", 
                     level = level, adjust = adjust, k = k)
    pred.lwr <- pred[, "lwr"]
    pred.upr <- pred[, "upr"]
    pred.ymin <- min(pred.lwr)
    pred.ymax <- max(pred.upr)
  }
  
  # Automatic limits for y-axis
  if (missing(ylim)) {
    if(interval == "prediction" || interval == "both") {
      ylim <- c(min(c(pred.ymin, yvals)), max(c(pred.ymax, yvals)))
    }
    if (interval == "confidence") {
      ylim <- c(min(c(conf.ymin, yvals)), max(c(conf.ymax, yvals)))
    } 
    if (interval == "none") {
      ylim <- c(min(c(fit.ymin, yvals)), max(c(fit.ymax, yvals)))
    }
  }
  
  # Plot data, mean response, etc.
  plot(xvals, yvals, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
       panel.first = if (hide) {
         if (shade) {
           # Draw (hidden) shaded prediction band
           if (interval == "prediction" || interval == "both") {
             polygon(c(xgrid[[1]], rev(xgrid[[1]])), c(pred.lwr, rev(pred.upr)), 
                     col = col.pred, border = border.pred, lty = lty.pred, 
                     lwd = lwd.pred)
           }
           # Draw (hidden) shaded confidence band
           if (interval == "confidence" || interval == "both") {
             polygon(c(xgrid[[1]], rev(xgrid[[1]])), c(conf.lwr, rev(conf.upr)), 
                     col = col.conf, border = border.conf, lty = lty.conf, 
                     lwd = lwd.conf)
           }
         } else {
           # Draw (hidden) unshaded prediction band
           if (interval == "prediction" || interval == "both") {
             lines(xgrid[[1]], pred.lwr, col = col.pred, lty = lty.pred, 
                   lwd = lwd.pred)
             lines(xgrid[[1]], pred.upr, col = col.pred, lty = lty.pred, 
                   lwd = lwd.pred)
           }
           # Draw (hidden) unshaded confidence band
           if (interval == "confidence" || interval == "both") {
             lines(xgrid[[1]], conf.lwr, col = col.conf, lty = lty.conf, 
                   lwd = lwd.conf)
             lines(xgrid[[1]], conf.upr, col = col.conf, lty = lty.conf, 
                   lwd = lwd.conf)
           }
         }
         # Draw (hidden) fitted response curve
         lines(xgrid[[1]], suppressWarnings(predict(object, newdata = xgrid)), 
               lty = lty.fit, lwd = lwd.fit, col = col.fit)
       } else NULL,
       panel.last = if(!hide) {
         if (shade) {
           # Draw shaded prediction band
           if (interval == "prediction" || interval == "both") {
             polygon(c(xgrid[[1]], rev(xgrid[[1]])), c(pred.lwr, rev(pred.upr)), 
                     col = col.pred, border = border.pred, lty = lty.pred, 
                     lwd = lwd.pred)
           }
           # Draw shaded confidence band
           if (interval == "confidence" || interval == "both") {
             polygon(c(xgrid[[1]], rev(xgrid[[1]])), c(conf.lwr, rev(conf.upr)), 
                     col = col.conf, border = border.conf, lty = lty.conf, 
                     lwd = lwd.conf)
           }
         } else {
           # Draw unshaded prediction band
           if (interval == "prediction" || interval == "both") {
             lines(xgrid[[1]], pred.lwr, col = col.pred, lty = lty.pred, 
                   lwd = lwd.pred)
             lines(xgrid[[1]], pred.upr, col = col.pred, lty = lty.pred, 
                   lwd = lwd.pred)
           }
           # Draw unshaded confidence band
           if (interval == "confidence" || interval == "both") {
             lines(xgrid[[1]], conf.lwr, col = col.conf, lty = lty.conf, 
                   lwd = lwd.conf)
             lines(xgrid[[1]], conf.upr, col = col.conf, lty = lty.conf, 
                   lwd = lwd.conf)
           }
         }
         # Draw fitted response curve
         lines(xgrid[[1]], suppressWarnings(predict(object, newdata = xgrid)), 
               lty = lty.fit, lwd = lwd.fit, col = col.fit)
       } else NULL, ...)
  
}


#' @rdname plotFit
#' @export
#' @method plotFit glm
plotFit.glm <- function(object, type = c("response", "link"),
                        interval = c("none", "confidence"), level = 0.95,
                        data, ..., shade = FALSE, extend.range = FALSE, 
                        hide = TRUE, 
                        col.conf = if (shade) grey(0.9) else "black", 
                        border.conf = col.conf, col.fit = "black", 
                        lty.conf = if (shade) 1 else 2, 
                        lty.fit = 1, lwd.conf = 1, lwd.fit = 1, n = 500, 
                        xlab, ylab, xlim, ylim) {
  
  # Extract data 
  if (!missing(data)) { 
    .data <- data 
  } else {
    .data <- eval(getCall(object)$data, envir = parent.frame())
  }
  if (is.null(.data)) {  # throw error if no data are found
    stop(paste("No data to plot. If the", class(object), "object does not",
               "contain a data component then a data frame containing the", 
               "variables in the model must be supplied through the 'data'", 
               "argument to plotFit."))
  }
  
  # Extract variable names and values
  xname <- intersect(all.vars(formula(object)[[3]]), colnames(.data)) 
  yname <- all.vars(formula(object)[[2]])
  if (length(xname) != 1) stop("Only one independent variable allowed.")
  xvals <- .data[, xname]
  
  # NOTE:
  #
  # For binomial and quasibinomial families the response can also be 
  # specified as a factor (when the first level denotes failure and all others 
  # success) or as a two-column matrix with the columns giving the numbers of 
  # successes and failures.
  if (family(object)$family %in% c("binomial", "quasibinomial")) {
    if (length(yname) == 1) {
      yvals <- with(.data, eval(formula(object)[[2]]))
    } else {
      ymat <- .data[, yname]
      yvals <- ymat[, 1] / ymat[, 2]
    }
  } else {
    if (length(yname) != 1) stop("Only one dependent variable allowed.")
    yvals <- with(.data, eval(formula(object)[[2]]))
  }
  
  # FIXME: Produce -Inf/Inf for binomial model for yvals = 0/1
  type <- match.arg(type)
  if (type == "link") yvals <- family(object)$linkfun(yvals)
  
  
  # Plot limits, labels, etc.
  if (missing(xlim)) xlim <- range(xvals)  # default limits for x-axis
  xgrid <- if (extend.range) {  # the x values at which to evaluate
    list(seq(from = extendrange(xlim)[1], to = extendrange(xlim)[2], 
             length = n))
  } else {
    list(seq(from = xlim[1], to = xlim[2], length = n))
  }
  names(xgrid) <- xname
  if (missing(xlab)) xlab <- xname  # default label for x-axis
  if (missing(ylab)) {
    ylab <- paste0("Prediction (", type, " scale)")  # default label for y-axis
  }
  
  # Confidence interval for mean response
  interval = match.arg(interval)
  if (interval == "confidence") {
    
    # Mean response and pointwise confidence intervals
    conf <- predict(object, newdata = xgrid, type = "link", se.fit = TRUE)
    lwr_link <- conf$fit - conf$se.fit * qnorm((level+1) / 2)
    upr_link <- conf$fit + conf$se.fit * qnorm((level+1) / 2)
    
    if (type == "response") {
      mean_resp <- family(object)$linkinv(conf$fit)
      lwr_conf <- family(object)$linkinv(lwr_link)
      upr_conf <- family(object)$linkinv(upr_link)
    } else {
      mean_resp <- conf$fit
      lwr_conf <- lwr_link
      upr_conf <- upr_link
    }

    # Determine default ylim for plot
    if (missing(ylim)) {
      ylim <- c(min(c(min(lwr_conf), yvals)), max(c(max(upr_conf), yvals)))
    }
    
  } else {
    
    # Mean response
    mean_resp <- unname(predict(object, newdata = xgrid, type = type))
    
    # Determine default ylim for plot
    if (missing(ylim)) {
      ylim <- c(min(c(min(mean_resp), yvals)), max(c(max(mean_resp), yvals)))
    }
    
  }
  
  # Plot data, mean response, etc.
  plot(xvals, yvals, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
       panel.first = if (hide) {
         if (shade) {
           # Draw (hidden) shaded confidence band
           if (interval == "confidence") {
             polygon(c(xgrid[[1]], rev(xgrid[[1]])), c(lwr_conf, rev(upr_conf)), 
                     col = col.conf, border = border.conf, lty = lty.conf, 
                     lwd = lwd.conf)
           }
         } else {
           # Draw (hidden) unshaded confidence band
           if (interval == "confidence") {
             lines(xgrid[[1]], lwr_conf, col = col.conf, lty = lty.conf, 
                   lwd = lwd.conf)
             lines(xgrid[[1]], upr_conf, col = col.conf, lty = lty.conf, 
                   lwd = lwd.conf)
           }
         }
         # Draw (hidden) fitted response curve
         lines(xgrid[[1]], mean_resp, lty = lty.fit, lwd = lwd.fit, 
               col = col.fit)
       } else NULL,
       panel.last = if (!hide) {
         if (shade) {
           # Draw shaded confidence band
           if (interval == "confidence") {
             polygon(c(xgrid[[1]], rev(xgrid[[1]])), c(lwr_conf, rev(upr_conf)), 
                     col = col.conf, border = border.conf, lty = lty.conf, 
                     lwd = lwd.conf)
           }
         } else {
           # Draw unshaded confidence band
           if (interval == "confidence") {
             lines(xgrid[[1]], lwr_conf, col = col.conf, lty = lty.conf, 
                   lwd = lwd.conf)
             lines(xgrid[[1]], upr_conf, col = col.conf, lty = lty.conf, 
                   lwd = lwd.conf)
           }
         }
         # Draw fitted response curve
         lines(xgrid[[1]], mean_resp, lty = lty.fit, lwd = lwd.fit, 
               col = col.fit)
       } else NULL, ...)
  
}


#' @rdname plotFit
#' @export
#' @method plotFit rlm
plotFit.rlm <- function(object, data, ..., extend.range = FALSE, hide = TRUE,
                        col.fit = "black", lty.fit = 1, lwd.fit = 1, n = 500, 
                        xlab, ylab, xlim, ylim) {
  
  # Extract data 
  if (!missing(data)) { 
    .data <- data 
  } else {
    .data <- eval(getCall(object)$data, envir = parent.frame())
  }
  if (is.null(.data)) {  # throw error if no data are found
    stop(paste("No data to plot. If the", class(object), "object does not",
               "contain a data component then a data frame containing the", 
               "variables in the model must be supplied through the 'data'", 
               "argument to plotFit."))
  }
  
  # Extract variable names and values
  xname <- intersect(all.vars(formula(object)[[3]]), colnames(.data)) 
  yname <- all.vars(formula(object)[[2]])
  if (length(xname) != 1) stop("Only one independent variable allowed.")
  if (length(yname) != 1) stop("Only one dependent variable allowed.")
  xvals <- .data[, xname]
  yvals <- with(.data, eval(formula(object)[[2]]))
  
  # Plot limits, labels, etc.
  if (missing(xlim)) xlim <- range(xvals)  # default limits for x-axis
  xgrid <- if (extend.range) {  # the x values at which to evaluate
    list(seq(from = extendrange(xlim)[1], to = extendrange(xlim)[2], 
             length = n))
  } else {
    list(seq(from = xlim[1], to = xlim[2], length = n))
  }
  names(xgrid) <- xname
  if (missing(xlab)) xlab <- xname  # default label for x-axis
  if (missing(ylab)) ylab <- yname  # default label for y-axis
  
  fitvals <- predict(object, newdata = xgrid)
  fit.ymin <- min(fitvals)
  fit.ymax <- max(fitvals)
  ylim <- c(min(c(fit.ymin, yvals)), max(c(fit.ymax, yvals)))

  
  # Plot data, mean response, etc.
  plot(xvals, yvals, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
       panel.first = if (hide) {  # draw points last
         # Draw (hidden) fitted response curve
         lines(xgrid[[1]], suppressWarnings(predict(object, newdata = xgrid)), 
               lty = lty.fit, lwd = lwd.fit, col = col.fit)  
       } else NULL,
       panel.last = if(!hide) {  # draw points first
         # Draw fitted response curve
         lines(xgrid[[1]], suppressWarnings(predict(object, newdata = xgrid)), 
               lty = lty.fit, lwd = lwd.fit, col = col.fit)  
       } else NULL, ...)
  
}

#' @rdname plotFit
#' @export
#' @method plotFit lqs
plotFit.lqs <- function(object, data, ..., extend.range = FALSE, hide = TRUE,
                        col.fit = "black", lty.fit = 1, lwd.fit = 1, n = 500, 
                        xlab, ylab, xlim, ylim) {
  
  # Extract data 
  if (!missing(data)) { 
    .data <- data 
  } else {
    .data <- eval(getCall(object)$data, envir = parent.frame())
  }
  if (is.null(.data)) {  # throw error if no data are found
    stop(paste("No data to plot. If the", class(object), "object does not",
               "contain a data component then a data frame containing the", 
               "variables in the model must be supplied through the 'data'", 
               "argument to plotFit."))
  }
  
  # Extract variable names and values
  xname <- intersect(all.vars(formula(object)[[3]]), colnames(.data)) 
  yname <- all.vars(formula(object)[[2]])
  if (length(xname) != 1) stop("Only one independent variable allowed.")
  if (length(yname) != 1) stop("Only one dependent variable allowed.")
  xvals <- .data[, xname]
  yvals <- with(.data, eval(formula(object)[[2]]))
  
  # Plot limits, labels, etc.
  if (missing(xlim)) xlim <- range(xvals)  # default limits for x-axis
  xgrid <- if (extend.range) {  # the x values at which to evaluate
    list(seq(from = extendrange(xlim)[1], to = extendrange(xlim)[2], 
             length = n))
  } else {
    list(seq(from = xlim[1], to = xlim[2], length = n))
  }
  names(xgrid) <- xname
  if (missing(xlab)) xlab <- xname  # default label for x-axis
  if (missing(ylab)) ylab <- yname  # default label for y-axis
  
  fitvals <- predict(object, newdata = xgrid)
  fit.ymin <- min(fitvals)
  fit.ymax <- max(fitvals)
  ylim <- c(min(c(fit.ymin, yvals)), max(c(fit.ymax, yvals)))
  
  
  # Plot data, mean response, etc.
  plot(xvals, yvals, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim,
       panel.first = if (hide) {  # draw points last
         # Draw (hidden) fitted response curve
         lines(xgrid[[1]], suppressWarnings(predict(object, newdata = xgrid)), 
               lty = lty.fit, lwd = lwd.fit, col = col.fit)  
       } else NULL,
       panel.last = if(!hide) {  # draw points first
         # Draw fitted response curve
         lines(xgrid[[1]], suppressWarnings(predict(object, newdata = xgrid)), 
               lty = lty.fit, lwd = lwd.fit, col = col.fit)  
       } else NULL, ...)
  
}

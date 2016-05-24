#' @name ciplot
#' @keywords barplot error bars
#' @author Sven E. Templer
#' @title Barplot with confindence intervals
#' @description 
#' Create barplots of a list of numeric values and error bars
#' according to the confidence interval, standard deviation,
#' interquartile range, etc.
#' @param x List of numeric values
#' @param ... Arguments forwarded to barplot in default method.
#' @param ylim A range for the y-axis limits.
#' @param height.fun Function to apply on each list object 
#' to calculate the height of the bars from.
#' @param height.args Arguments forwarded to height.fun, as a named list.
#' @param error.fun Function to calculate the error size. See also details.
#' @param error.args Arguments forwarded to error.fun, as a named list.
#' @param arrows.args Arguments forwarded to arrows, as a named list.
#' @param na.rm Logical, remove missing values.
#' @details 
#' Example for quantiles:\cr
#' \code{interquartile <- function(x) \{quartile(x,.75)-mean(x)\}}\cr
#' \code{quantileQ <- function(x, q) \{abs(quartile(x,q[1])-mean(x))\}}

#' @rdname ciplot
#' @export
ciplot <- function (x, ...) { UseMethod("ciplot") }

#' @rdname ciplot
#' @export
ciplot.default <- function (x, ..., ylim,
                            height.fun = mean, height.args = list(), 
                            error.fun = confint, error.args = list(),
                            arrows.args = list(code = 3, angle = 90), 
                            na.rm = TRUE) {
  
  if (na.rm)
    x <- lapply(x, function (y) y[!is.na(y)])
  h <- sapply(x, function(y) do.call(height.fun, c(list(y), height.args)))
  e <- sapply(x, function(y) do.call(error.fun, c(list(y), error.args)))
  e.min <- h - e
  e.max <- h + e
  if (missing(ylim)) {
    al <- min(-1e-8, e.min * 1.1)
    ar <- max(1e-8, e.max * 1.1)
    ylim <- c(al, ar)
  }
    
  b <- barplot(h, ylim = ylim, ...)
  b <- sapply(b, as.vector)
  arrows.args <- c(list(x0 = b, x1 = b, y0 = e.min, y1 = e.max), 
                   arrows.args)
  null <- do.call(arrows, arrows.args)
  invisible(list(height=h, error=e))
}


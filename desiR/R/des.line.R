#' @title Plots a desirability function on an existing graph
#'
#' @description Plots any of the desirability functions on top of a graph,
#' usually a histogram or density plot.
#'
#' @details This function can be used to visualise how the desirabilities are
#' mapped from the raw data to a 0-1 scale, which can help select suitable cut
#' points. The scale of the y-axis has a minimum of 0 and a maximum of 1.
#'
#' WARNING: If you set xlim values for the histogram or density plot, then you must pass the same xlim values to des.line; otherwise the data and desirability function (plotted line) will be misaligned. If xlim is not set, then the same default values will be used for the data and the function. 
#'
#' @param x Vector of numeric or integer values.
#' @param des.func Name of the desirability function to plot (in quotes).
#' @param des.args A vector of named arguments for the chosen desirability
#' function.
#' @param ... Arguments for the plotting function (e.g. xlim, lwd, lty).
#'
#' @return Plotted values of the desirability function.
#' @seealso \code{\link{d.low}}, \code{\link{d.high}}, \code{\link{d.central}},
#' \code{\link{d.ends}}, \code{\link{d.4pl}}
#' 
#' @examples
#' set.seed(1)
#' x1 <- rnorm(100, 10, 2)
#' hist(x1, breaks=10)
#' des.line(x1, "d.high", des.args=c(cut1=10, cut2=11))
#' des.line(x1, "d.high", des.args=c(cut1=10, cut2=11,
#' des.min=0.1, scale=0.5))



des.line <- function(x, des.func, des.args, ...){

  # check function selected correctly (one must be specified)
  if (!des.func %in% c("d.low","d.high","d.central", "d.ends",
                       "d.4pl", "d.rank")) {
    stop("\ndes.func must be one of 'd.high', 'd.low',
         'd.central', 'd.ends', 'd.4pl', or 'd.rank'\n")
  }

  # make values for plotting
  vals <- seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE), length.out=300)

  # add lines to existing graph
  graphics::par(new=TRUE)
  
  if (des.func=="d.high") {
   graphics:: plot(do.call(d.high, c(list(x=vals), des.args) ) ~ vals,
         ylim=c(0,1), type="l",  yaxt="n", bty="n", xaxt="n",
         ylab="", xlab="", ...)
  }

  
  if (des.func=="d.low") {
    graphics::plot(do.call(d.low, c(list(x=vals), des.args) ) ~ vals,
         ylim=c(0,1), type="l",  yaxt="n", bty="n", xaxt="n",
         ylab="", xlab="", ...)
  }

  
  if (des.func=="d.central") {
    graphics::plot(do.call(d.central, c(list(x=vals), des.args) ) ~ vals,
         ylim=c(0,1), type="l",  yaxt="n", bty="n", xaxt="n",
         ylab="", xlab="", ...)
  }

  
  if (des.func=="d.ends") {
    graphics::plot(do.call(d.ends, c(list(x=vals), des.args) ) ~ vals,
         ylim=c(0,1), type="l",  yaxt="n", bty="n", xaxt="n",
         ylab="", xlab="", ...)
  }

    
  if (des.func=="d.4pl") {
    graphics::plot(do.call(d.4pl, c(list(x=vals), des.args) ) ~ vals,
         ylim=c(0,1), type="l",  yaxt="n", bty="n", xaxt="n",
         ylab="", xlab="", ...)
  }

  if (des.func=="d.rank") {
      graphics::plot(do.call(d.rank, c(list(x=vals), des.args) ) ~ vals,
         ylim=c(0,1), type="l",  yaxt="n", bty="n", xaxt="n",
         ylab="", xlab="", ...)
  }

}

#' Plotting methods
#' 
#' Use this function to produce plots for an object of class \code{sim_setup}.
#' 
#' @param x a \code{sim_setup}
#' @param y will be ignored
#' @param ... Arguments to be passed to \code{\link[graphics]{plot}}.
#' 
#' @export
#' @method plot sim_setup
#' @rdname plot.sim_setup
#' @seealso \code{\link[saeSim]{autoplot}}
plot.sim_setup <- function(x, y, ...) {
  # Get some data
  dat <- as.data.frame(x)
  # get name for y-Axis; default to "y" 
  yAxis <- if(all(names(dat) != "y")) 
    names(dat)[!grepl("id", names(dat))][1] else "y"
  # get name for x-axis
  xAxis <- names(dat)[!grepl(paste("id", yAxis, sep = "|"), names(dat))][1]
  
  plot(y=dat[[yAxis]], x = dat[[xAxis]], xlab = xAxis, ylab = yAxis, ...)
}
#' Autoplot method
#' 
#' Use this function to produce plots for an object of class \code{sim_setup} and you like to have plots based on ggplot2. At this time it is a ggplot2 implementation which mimics the behavior of \code{\link{smoothScatter}} without all the options.
#' 
#' @param object a sim_setup
#' @param x character of variable name in the data on the x-axis
#' @param y character of variable name in the data on the y-axis
#' @param ... is not used
#' 
#' @export
#' @method autoplot sim_setup
#' @rdname autoplot
#' @examples
#' \dontrun{
#' autoplot(sim_base_lm())
#' }
autoplot.sim_setup <-  function(object, x = "x", y = "y", ...) {
  # Get some data
  dat <- as.data.frame(object)
  # get name for y-Axis; default to "y" 
  yAxis <- if(all(names(dat) != y)) {
    warning("Variable for y axis not found, picking one automatically.")
    names(dat)[!grepl("id", names(dat))][1]
  } else y
  # get name for x-axis
  xAxis <- if(all(names(dat) != x)) {
    warning("Variable for x axis not found, picking one automatically.")
    names(dat)[!grepl(paste("id", yAxis, sep = "|"), names(dat))][1]
    } else x
  
  # Make the plot
  ggplot(dat, aes_string(x = xAxis, y = yAxis)) + 
    stat_density2d(geom="tile", aes_string(fill="..density..^0.25", alpha="1"), contour=FALSE) + 
    geom_point(alpha = 0.1, size = 0.5) +
    stat_density2d(geom="tile", 
                   aes_string(fill="..density..^0.25", 
                       alpha="ifelse(..density..^0.25<0.4,0,1)"), 
                   contour=FALSE) +
    scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256)) + 
    theme_classic() + theme(legend.position = "none")
}

#' @name autoplot
#' @export autoplot
#' @rdname autoplot
NULL
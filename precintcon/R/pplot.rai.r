#' @include precintcon.plot.rai.r
NULL

#' @name pplot.rai
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com} 
#' @aliases precintcon.plot.rai pplot.rai 
#' @title Plot Rainfall Anomaly Index 
#' @description Plots the Rainfall Anomaly Index of a precipitation serie. 
#' @usage pplot.rai(\dots, granularity = "m", xlab = "Month", ylab = "RAI", 
#' 					 ylim = c(-3,3), legend = NULL, fontsize = 10, 
#'                 axis.text.color = "black", export = FALSE, 
#' 					 export.name = "rai_plot.png", width = 8.6, 
#'					    height = 7.5, units = "cm") 
#' @param \dots a set of daily or monthly precipitation series.
#' @param granularity the granularity applied for calculating the rainfall 
#' anomaly index, which should be either "a" for annual granularity or "m" for 
#' monthly granularity ". (Default value: "m")
#' @param xlab the text for the x axis. (Default value: "Month")
#' @param ylab the text for the y axis. (Default value: "RAI")
#' @param ylim the limits of the y axis. (Default value: c(-3, 3))
#' @param legend the text vector for the legend items. If NULL the legends will 
#' be equals to the variable names. (Default value: NULL)
#' @param fontsize the font size value in pt. (Default value: 10)
#' @param axis.text.color the legend colors. (Default value: "black")
#' @param export the logical value for defining whether the graph should be export 
#' to a file or not. (Default value: FALSE)
#' @param export.name the text for defining the exported file name. It is only used 
#' if export = TRUE. (Default value: "rai_plot.png")
#' @param width the number for defining the exported graph width. It is only used if 
#' export = TRUE. (Default value: 8.6)
#' @param height the number for defining the exported graph height. It is only used 
#' if export = TRUE. (Default value: 7.5)
#' @param units the text for defining the units of the height and width parameters. 
#' It is only used if export = TRUE. (Default value: "cm")
#' @seealso \code{\link{read.data}}
#' @examples 
#' ##
#' # Loading the daily precipitation serie.
#' data(daily)
#' 
#' ##
#' # Performing the a set of statistical analysis
#' pplot.rai(daily, granularity = "m")
#' @references Rooy, M. P. van. A Rainfall anomaly index independent of time and space, Notos. v.14, p.43-48, 1965.
#' @keywords rainfall anomaly precipitation 
#' @export
pplot.rai <- function(
  ..., 
	granularity     = "m",
	xlab            = "Month",
	ylab            = "RAI", 
	ylim            = c(-3,3),
	legend          = NULL,
	fontsize        = 10, 
	axis.text.color = "black", 
	export          = FALSE, 
	export.name     = "rai_plot.png", 
	width           = 8.6, 
  height          = 7.5, 
  units            = "cm"
) {
   precintcon.plot.rai(..., granularity = granularity, xlab = xlab, ylab = ylab, 
         ylim = ylim, legend = legend, fontsize = fontsize, 
         axis.text.color = axis.text.color, export = export, 
         export.name = export.name, width = width, height = height, units = units,
         args = as.character(match.call()[1:length(list(...))+1]))
}
#' @include precintcon.plot.histogram.r
NULL

#' @name pplot.histogram
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com} 
#' @aliases precintcon.plot.histogram pplot.histogram 
#' @title Plot histogram
#' @description Plots the histogram of a precipitation serie. 
#' @usage pplot.histogram(\dots, density = FALSE, xlab = "Precipitation (mm)", 
#' ylab = "Frequency", legend.title = "Legend", 
#' legend = NULL, fontsize = 10, axis.text.color = "black", 
#' export = FALSE, export.name = "histogram_plot.png", 
#' width = 10, height = 10, units = "cm") 
#' @param \dots a set of daily or monthly precipitation series.
#' @param density the logical value for defining whether the graph should be 
#' plotted with bars or lines. (Default value: FALSE)
#' @param xlab the text for the x axis. (Default value: "Precipitation (mm)")
#' @param ylab the text for the y axis. (Default value: "Frequency")
#' @param legend.title the text for the legend title. (Default value: "Legend")
#' @param legend the text vector for the legend items. If NULL the legends will 
#' be equals to the variable names. (Default value: NULL)
#' @param fontsize the font size value in pt. (Default value: 10)
#' @param axis.text.color the legend colors. (Default value: "black")
#' @param export the logical value for defining whether the graph should be 
#' export to a file or not. (Default value: FALSE)
#' @param export.name the text for defining the exported file name. It is only 
#' used if export = TRUE. (Default value: "histogram_plot.png")
#' @param width the number for defining the exported graph width. It is only 
#' used if export = TRUE. (Default value: 10)
#' @param height the number for defining the exported graph height. It is only 
#' used if export = TRUE. (Default value: 10)
#' @param units the text for defining the units of the height and width 
#' parameters. It is only used if export = TRUE. (Default value: "cm")
#' @seealso
#' \code{\link{read.data}}
#' @examples 
#' ##
#' # Loading the daily precipitation serie.
#' data(daily)
#' 
#' ##
#' # Performing the a set of statistical analysis
#' pplot.histogram(daily)
#' @keywords histogram precipitation
#' @export 
pplot.histogram <- function(
   ..., 
	density         = FALSE, 
	xlab            = "Precipitation (mm)",
	ylab            = "Frequency",
  legend.title    = "Legend", 
  legend          = NULL,
	fontsize        = 10,
	axis.text.color = "black", 
	export          = FALSE, 
	export.name     = "histogram_plot.png", 
	width           = 10, 
	height          = 10, 
	units           = "cm"
) {
   precintcon.plot.histogram(..., density = density,
         xlab = xlab, ylab = ylab, legend.title = legend.title, legend = legend, 
         fontsize = fontsize, axis.text.color = axis.text.color, export = export, 
         export.name = export.name, width = width, height = height, units = units,
         args = as.character(match.call()[1:length(list(...))+1]))	
}

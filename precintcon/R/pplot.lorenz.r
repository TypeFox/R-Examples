#' @include precintcon.plot.lorenz.r
NULL

#' @name pplot.lorenz
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com}
#' @aliases precintcon.plot.lorenz pplot.lorenz
#' @title Plot Lorenz's curve
#' @description Plots the Lorenz's curve of a precipitation serie.
#' @usage pplot.lorenz(\dots, interval = 1, grouped = FALSE, 
#' xlab = expression(sum(n[i]), i==1), ylab = expression(sum(P[i]), i==1), 
#' legend.title = "Legend", legend = NULL, fontsize = 10, 
#' axis.text.color = "black", export = FALSE, export.name = "lorenz_plot.png", 
#' width = 8.6, height = 7.5, units = "cm")
#' @param \dots a set of daily precipitation series.
#' @param interval the interval in millimeters applied for calculating the 
#' Lorenz's curve. (Default value: 1)
#' @param grouped the logical value for defining whether all series should 
#' be plotted in the same graph or not. (Default value: FALSE)
#' @param xlab the text for the x axis. 
#' (Default value: expression(sum(n[i]), i==1)
#' @param ylab the text for the y axis. 
#' (Default value: expression(sum(P[i]), i==1))
#' @param legend.title the text for the legend title. 
#' (Default value: "Legend")
#' @param legend the text vector for the legend items. 
#' If NULL the legends will be equals to the variable names. (Default value: NULL)
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
#' \code{\link{ci}} 
#' \code{\link{read.data}}
#' \code{\link{as.daily}}
#' @examples
#' ##
#' # Loading the daily precipitation serie.
#' data(daily)
#'
#' ##
#' # Performing the a set of statistical analysis
#' pplot.lorenz(daily, interval = 1)
#' @keywords lorenz's curve precipitation
#' @export
pplot.lorenz <- function(
  ..., 
  interval        = 1,
	grouped         = FALSE, 
	xlab            = expression(sum(n[i]), i==1),
	ylab            = expression(sum(P[i]), i==1), 
	legend.title    = "Legend",
	legend          = NULL,
	fontsize        = 10, 
	axis.text.color = "black", 
	export          = FALSE, 
	export.name     = "lorenz_plot.png", 
  width           = 8.6, 
  height          = 7.5, 
  units           = "cm"
)
precintcon.plot.lorenz(..., 
  interval        = interval, 
  grouped         = grouped, 
  xlab            = xlab, 
  ylab            = ylab, 
  legend.title    = legend.title, 
  legend          = legend, 
  fontsize        = fontsize, 
  axis.text.color = axis.text.color, 
  export          = export, 
  export.name     = export.name, 
  width           = width, 
  height          = height, 
  units           = units,
  args            = as.character(match.call()[1:length(list(...))+1]))
#' @include precintcon.plot.deciles.r
NULL

#' @name pplot.deciles
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com}
#' @aliases precintcon.plot.deciles pplot.deciles
#' @title Plot deciles
#' @description Plots the deciles of a precipitation serie.
#' @usage pplot.deciles(\dots, ylab = "Precipitation", 
#' legend.title = "Legend", legend = NULL, 
#' fontsize = 10, axis.text.color = "black", 
#' export = FALSE, export.name = "deciles_plot.png",
#' width = 8.6, height = 7.5, units = "cm", grouped = FALSE)
#' @param \dots a set of daily or monthly precipitation serie.
#' @param ylab the text for the y axis. (Default value: "Precipitation")
#' @param legend.title the text for the legend title. (Default value: "Legend")
#' @param legend the text vector for the legend items. If NULL the legends will 
#' be equals to the variable names. (Default value: NULL)
#' @param fontsize the font size value in pt. (Default value: 10)
#' @param axis.text.color the legend colors. (Default value: "black")
#' @param export the logical value for defining whether the graph should be 
#' export to a file or not. (Default value: FALSE)
#' @param export.name the text for defining the exported file name. It is only 
#' used if export = TRUE. (Default value: "deciles_plot.png")
#' @param width the number for defining the exported graph width. It is only used 
#' if export = TRUE. (Default value: 8.6)
#' @param height the number for defining the exported graph height. It is only 
#' used if export = TRUE. (Default value: 7.5)
#' @param units the text for defining the units of the height and width 
#' parameters. It is only used if export = TRUE. (Default value: "cm")
#' @param grouped the logical value for defining whether the graphs should be plotted in
#' group. 
#' @seealso 
#' \code{\link{deciles}}
#' \code{\link{read.data}}
#' @examples
#' ## Loading the monthly precipitation serie.
#' #
#' data(monthly)
#' 
#' ## Performing the a set of statistical analysis
#' #
#' pplot.deciles(monthly)
#' @keywords deciles precipitation
#' @export 
pplot.deciles <- function(...,
	ylab            = "Precipitation",
	legend.title    = "Legend",
	legend          = NULL,
	fontsize        = 10, 
	axis.text.color = "black", 
	export          = FALSE, 
	export.name     = "deciles_plot.png", 
	width           = 8.6, 
	height          = 7.5, 
	units           = "cm",
	grouped         = FALSE
)
	precintcon.plot.deciles(..., 
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
  grouped         = grouped,
	args            = as.character(match.call()[1:length(list(...))+1]))
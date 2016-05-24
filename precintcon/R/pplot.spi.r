#' @include precintcon.plot.spi.r
NULL

#' @name pplot.spi
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com} 
#' @aliases precintcon.plot.spi pplot.spi 
#' @title Plot Standardized Precipitation Index 
#' @description Plots the Standardized Precipitation Index of a precipitation 
#' serie. 
#' @usage pplot.spi(\dots, period = 3, distribution = "Gamma", xlab = "Months",
#' ylab = "SPI", ylim = c(-3,3), legend = NULL, fontsize = 10, 
#' axis.text.color = "black", export = FALSE, export.name = "spi_plot.png", 
#' width = 8.6, height = 7.5, units = "cm") 
#' @param \dots a set of daily or monthly precipitation series.
#' @param period the number of months to be aggregated for calculating 
#'                 the standardized precipitation index. (Default value: 3)
#' @param distribution it has no effect yet. (Default value: "Gamma")
#' @param xlab the text for the x axis. (Default value: "Months")
#' @param ylab the text for the y axis. (Default value: "SPI")
#' @param ylim the limits of the y axis. (Default value: c(-3, 3))
#' @param legend the text vector for the legend items. If NULL the legends will 
#' be equals to the variable names. (Default value: NULL)
#' @param fontsize the font size value in pt. (Default value: 10)
#' @param axis.text.color the legend colors. (Default value: "black")
#' @param export the logical value for defining whether the graph should be export 
#' to a file or not. (Default value: FALSE)
#' @param export.name the text for defining the exported file name. It is only used 
#' if export = TRUE. (Default value: "spi_plot.png")
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
#' pplot.spi(daily)
#' @keywords standardized precipitation index precipitation
#' @export 
pplot.spi <- function(..., 
	period          = 3, 
  distribution    = "Gamma",
	xlab            = "Months",
	ylab            = "SPI", 
	ylim            = c(-3,3),
	legend          = NULL,
	fontsize        = 10, 
	axis.text.color = "black", 
	export          = FALSE, 
	export.name     = "spi_plot.png", 
	width           = 8.6, 
  height          = 7.5, 
  units           = "cm"
)

precintcon.plot.spi(..., 
  period          = period, 
  distribution    = distribution, 
  xlab            = xlab, 
  ylab            = ylab,
  ylim            = ylim, 
  legend          = legend, 
  fontsize        = fontsize, 
  axis.text.color = axis.text.color, 
  export          = export, 
  export.name     = export.name, 
  width           = width, 
  height          = height, 
  units           = units,
  args            = as.character(match.call()[1:length(list(...))+1]))
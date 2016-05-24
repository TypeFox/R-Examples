#' @include precintcon.plot.pci.r
NULL

#' @name pplot.pci
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com}
#' @aliases precintcon.plot.pci pplot.pci
#' @title Plot Precipitation Concentration Index
#' @description Plots the Precipitation Concentration Index of a precipitation 
#' serie.
#' @usage pplot.pci(\dots, xlab = "Years", ylab = "PCI", legend = NULL, 
#' fontsize = 10, axis.text.color = "black", export = FALSE, 
#' export.name = "pci_plot.png", width = 10, height = 10, units = "cm")
#' @param \dots a set of daily precipitation series.
#' @param xlab the text for the x axis. (Default value: "Years"
#' @param ylab the text for the y axis. (Default value: "PCI")
#' @param legend the text vector for the legend items. If NULL the legends 
#' will be equals to the variable names. (Default value: NULL)
#' @param fontsize the font size value in pt. (Default value: 10)
#' @param axis.text.color the legend colors. (Default value: "black") 
#' @param export the logical value for defining whether the graph should be 
#' export to a file or not. (Default value: FALSE)
#' @param export.name the text for defining the exported file name. It is only 
#' used if export = TRUE. (Default value: "pci_plot.png") 
#' @param width the number for defining the exported graph width. It is only 
#' used if export = TRUE. (Default value: 10)
#' @param height the number for defining the exported graph height. It is only 
#' used if export = TRUE. (Default value: 10) 
#' @param units the text for defining the units of the height and width 
#' parameters. It is only used if export = TRUE. (Default value: "cm")
#' @seealso \code{\link{read.data}}
#' @examples
#' ##
#' # Loading the daily precipitation serie.
#' data(daily)
#' 
#' ##
#' # Performing the a set of statistical analysis
#' pplot.pci(daily)
#' @keywords precipitation concentration index
#' @export
pplot.pci <- function(
  ..., 
	xlab            = "Years", 
  ylab            = "PCI", 
  legend          = NULL,
  fontsize        = 10, 
	axis.text.color = "black",  
	export          = FALSE, 
  export.name     = "pci_plot.png", 
	width           = 10, 
  height          = 10, 
  units           = "cm"
)
precintcon.plot.pci(..., 
  xlab            = xlab, 
  ylab            = ylab, 
  legend          = legend, 
  fontsize        = fontsize, 
  axis.text.color = axis.text.color, 
  export          = export, 
  export.name     = export.name, 
  width           = width, 
  height          = height, 
  units           = units,
  args            = as.character(match.call()[1:length(list(...))+1]))
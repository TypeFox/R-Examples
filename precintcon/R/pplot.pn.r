#' @include precintcon.plot.pn.r
NULL

#' @name pplot.pn
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com}
#' @aliases precintcon.plot.pn pplot.pn
#' @title Plot Percent of Normal
#' @description Plots the Percent of Normal of a precipitation serie.
#' @usage pplot.pn(\dots, interval = 30, scale = "a", xlab = NA, ylab = "PN", 
#'          fontsize = 10, axis.text.color = "black", legend = NULL, 
#'          export = FALSE, export.name = "pn_plot.png", width = 10, 
#'          height = 10, units = "cm")
#' @param \dots a set of daily or monthly precipitation serie.
#' @param interval the number of months applied for calculating the percentage 
#' of normal.
#' @param scale the scale used for calculating the percentage of normal, 
#'       which should be either "w" for weak (not supported yet), 
#'       "m" for month, "s" for season, or "d" for decades.
#' @param xlab the text for the x axis. (Default value: NA)
#' @param ylab the text for the y axis. (Default value: "PN")
#' @param legend the text vector for the legend items. If NULL the legends will 
#' be equals to the variable names. (Default value: NULL)
#' @param fontsize the font size value in pt. (Default value: 10)
#' @param axis.text.color the legend colors. (Default value: "black")
#' @param export the logical value for defining whether the graph should be export 
#' to a file or not. (Default value: FALSE)
#' @param export.name the text for defining the exported file name. It is only used 
#' if export = TRUE. (Default value: "pci_plot.png") 
#' @param width the number for defining the exported graph width. It is only used if 
#' export = TRUE. (Default value: 10)
#' @param height the number for defining the exported graph height. It is only used 
#' if export = TRUE. (Default value: 10)
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
#' pplot.pn(daily)
#' @keywords percent of normal precipitation
#' @export
pplot.pn <- function(..., 
  interval        = 30, 
  scale           = "a",
	xlab            = NA, 
  ylab            = "PN", 
  fontsize        = 10, 
	axis.text.color = "black", 
  legend          = NULL, 
	export          = FALSE, 
  export.name     = "pn_plot.png", 
	width           = 10, 
  height          = 10, 
  units            = "cm"
)

precintcon.plot.pn(..., 
  interval        = interval, 
  scale           = scale, 
  xlab            = xlab, 
  ylab            = ylab, 
  fontsize        = fontsize, 
  axis.text.color = axis.text.color, 
  legend          = legend, 
  export          = export, 
  export.name     = export.name, 
  width           = width, 
  height          = height, 
  units           = units,
  args            = as.character(match.call()[1:length(list(...))+1]))

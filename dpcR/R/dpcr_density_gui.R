#' Digital PCR Density Graphical User Interface
#' 
#' Launches graphical user interface that allows calculating density of
#' positive partitions distribution.
#' 
#' @section Warning : Any ad-blocking software may cause malfunctions.
#' @author Michal Burdukiewicz, Stefan Roediger.
#' @seealso \code{\link{dpcr_density}}.
#' @keywords hplot
#' @export dpcr_density_gui
dpcr_density_gui <- function()
  runApp(system.file("dpcr_density_gui", package = "dpcR"))

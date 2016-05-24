#' Digital PCR Report Graphical User Interface
#' 
#' Launches graphical user interface that generates reports from digital PCR 
#' data.
#' 
#' @section Warning : Any ad-blocking software may cause malfunctions.
#' @author Michal Burdukiewicz, Stefan Roediger.
#' @keywords hplot
#' @export dpcReport
dpcReport <- function()
  runApp(system.file("dpcReport", package = "dpcR"))

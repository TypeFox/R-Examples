#' Compare digital PCR runs - interactive presentation
#' 
#' Launches graphical user interface allowing multiple comparisons of simulated digital 
#' PCR reactions. 
#' 
#' @section Warning : Any ad-blocking software may be cause of malfunctions.
#' @author Michal Burdukiewicz, Stefan Roediger.
#' @seealso \code{\link{test_counts}}.
#' @keywords hplot
#' @export
test_counts_gui <- function()
  runApp(system.file("test_counts_gui", package = "dpcR"))
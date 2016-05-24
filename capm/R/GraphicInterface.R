#' Graphic interface to use some capm functions
#' @param set.func string to select a graphic interface for a set of functions to achieve a specific task. \code{SelectSamplingUnits} creates a graphic interface to select sampling units for pilot or final survey designs. \code{CalculateSampleSize} creates a graphic interface to calculate sample size and composition. The graphic interface created by \code{SurveyAnalysis} support functionality to analyise survey data. \code{SolveIASA} creates an interface to simulate population dynamics and to assess population-based interventions.
#' @return a graphic interface in a browser.
#' @details The graphic interfaces are created with the \code{shiny} package and thus will open in a browser.
#' @references \url{http://oswaldosantos.github.io/capm}
#' @export
#' 
#' @examples
#' # Uncomment the following line to open the graphic interface in a browser:
#' # GraphicInterface(set.func = 'SelectSamplingUnits')
#' 
GraphicInterface <- function(set.func) {
  if (Sys.info()['sysname'] == 'Windows') {
    set.func <- paste0('/', set.func)
  }
  runApp(paste0(system.file('shinyApps/', package = 'capm'), '/', set.func))
}
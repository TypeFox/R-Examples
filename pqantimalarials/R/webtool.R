#' @title
#' Start the pqantimalarials interactive web tool.
#' @description
#' This function starts the interactive web tool which
#' runs locally on the user's machine. If the browser doesn't
#' open automatically, the function provides the url so you can
#' manually open your browser and view the site. This function
#' normally does not return; interrupt R to stop the application
#' (usually by pressing Esc or Ctrl+C).
#' @seealso \code{\link{runApp}} which this function wraps.
#' @details
#' Within the web tool, the left panel allows users to set input
#' parameters that are then used to estimate the number of under-five
#' child deaths associated with poor-quality antimalarials in
#' sub-Saharan Africa according to the model presented in
#' Renschler et al. (2014). Users can browse through the site's tabs to
#' explore the outputs from uncertainty and sensitivity analyses that
#' were performed using their input settings. Users can download the
#' output data and visualizations as CSVs and PDFs respectively.
#' @examples
#' \dontrun{
#' # Start the web tool
#' webtool()
#' }
#' @references
#' J. Patrick Renschler, Kelsey Walters, Paul Newton,
#'   Ramanan Laxminarayan. "Estimated under-five deaths
#'   associated with poor-quality antimalarials in
#'   sub-Sarahan Africa". 2014.
#' @export

webtool <- function(){

  shiny::runApp(system.file('webtool', package = 'pqantimalarials'))
  
}
#' @details 
#' 
#' A light R wrapper around the elex command
#'   
#' elexr provides an R interface for the elex package 
#' (http://github.com/newsdev/elex/) to load election results from the 
#' Associated Press Elections API.
#'   
#' Rather than reinventing the wheel, this package simply runs \code{elex} in
#' a shell and loads the CSV output into a data frame.  Because of this,
#' you'll need Python and to install elex by following the instructions at 
#' \url{http://elex.readthedocs.org/en/latest/install.html}.
#'   
#' You'll also need an API key for the API Elections API.  You'll need to set
#' the \code{AP_API_KEY} environtment variable in your R session to be able to
#' retreive results.
#' 
#' @examples
#' library(elexr)
#' # Set your actual AP Elections API Key here. The string below is an example
#' # only and won't work.
#' Sys.setenv(AP_API_KEY = "AZA0AZ0aZ0ZA0az0AzAZ0AzazAzaZaZAZ")
#' tryCatch({
#'     ia_results <- results("02-01-2016")
#' },
#' error = function(err) {
#'     print(err)
#' })
#'   
"_PACKAGE"

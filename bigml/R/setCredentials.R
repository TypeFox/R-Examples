#' Set BigML API authentication credentials
#' @export 
#' @param username use the given username for all subsequent API requests
#' @param api_key use the given api key for all subsequent API requests
#' @template author
#' @details This function sets default username and api_key information for 
#'	subsequent BigML API access calls.  The relevent username and key are 
#'	stored in the R system environment variables.  So, it's also possible
#'	to set these variables by setting BIGMLUSER and BIGMLAPIKEY in an 
#'	.Renviron file.
#' @examples 
#' \dontrun{
#' # replace with your valid credentials:
#'	setCredentials('some_username', 'some_key')
#' }
setCredentials <-
function (username, api_key) 
{
    Sys.setenv(BIGML_USERNAME=username)
    Sys.setenv(BIGML_API_KEY=api_key)
}

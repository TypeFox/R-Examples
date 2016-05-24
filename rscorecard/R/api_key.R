#' Store Data.gov API key in system environment.
#'
#' This function stores your data.gov API key in the system environment
#' so that you only have to load it once at the start of the session.
#' If you set your key using \code{sc_key}, then you may omit
#' \code{api_key} parameter in the \code{sc_get} function.
#'
#' @param api_key Personal API key requested from
#'     \url{https://api.data.gov/signup} stored in a string.
#'
#' @examples
#' \dontrun{
#' sc_key('<API KEY IN STRING>')
#' }
#'
#' @section Obtain a key:
#' To obtain an API key, visit \url{https://api.data.gov/signup}.

#' @export
sc_key <- function(api_key) {

    ## check key
    if (nchar(api_key) != 40 || !is.character(api_key)) {
        stop('API key must be character string of length 40', call. = FALSE)
    }

    Sys.setenv(DATAGOV_API_KEY = api_key)
    message('DATAGOV_API_KEY environment variable now set. ' %+%
            'You may now use sc_get() without specifying a key.')
}




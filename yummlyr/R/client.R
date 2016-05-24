#' Save API credentials for later use
#'
#' This functions caches the credentials to avoid need for entering it when calling other functions
#' @param app_id application ID
#' @param app_key application key
#' @examples
#' # since not checking is preformed not to waste API calls
#' # it falls on the user to save correct information
#' save_yummly_credentials("APP_ID", "APP_KEY")
#' @export
save_yummly_credentials <- function(app_id, app_key) {
    if (app_id != "" && app_key != "") {
        assign("APP_ID", app_id, envir=auth_cache)
        assign("APP_KEY", app_key, envir=auth_cache)
    }
}

#' Get recipe from Yummly.com
#' 
#' This call is equivalent in functionality to a Yummly recipe page.
#' @param recipe_id recipe ID
#' @param app_id application ID
#' @param app_key application KEY
#' @note This function resembles viewing a recipe on Yummly.com
#' @references \itemize{
#'   \item Yummly Developer Guide \url{https://developer.yummly.com/documentation}
#' }
#' @examples
#' \dontrun{
#' # to request the response for French Onion Soup by Ree Drummond The Pioneer Woman 
#' # with id French-Onion-Soup-The-Pioneer-Woman-Cooks-_-Ree-Drummond-41364
#' get_recipe("French-Onion-Soup-The-Pioneer-Woman-Cooks-_-Ree-Drummond-41364")
#' }
#' @export
get_recipe <- function(recipe_id, app_id = auth_cache$APP_ID, app_key = auth_cache$APP_KEY) {
    if (is.null(app_id) || is.null(app_key)) {
        stop("APP_ID or APP_KEY is not set. Use save_yummly_credentials or supply appropriate arguments")
    }
    query <- sprintf("%s/%s?_app_id=%s&_app_key=%s", URL_GET,
                     recipe_id, app_id, app_key)
    content <- perform_query(utils::URLencode(query))
    jsonlite::fromJSON(content)
}
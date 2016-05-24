#' Fetch all available typeforms
#'
#' This function returns a data frame containing your typeforms and their
#' associated UIDs.
#' @importFrom jsonlite fromJSON
#' @param api Default \code{NULL}. Your private api key. If \code{api} is \code{NULL}
#' we use the environment variable \code{Sys.getenv("typeform_api")}.
#' @return A two column data frame.
#' @export
#' @examples
#' \dontrun{
#' api = "XXXXX"
#' get_all_typeforms(api)
#' }
get_all_typeforms = function(api=NULL) {
  api = get_api(api)
  url = paste0("https://api.typeform.com/v1/forms?key=", api)
  jsonlite::fromJSON(url)
}



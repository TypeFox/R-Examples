#' Extract Form Id From Form URL
#'
#' Extract Form Id From Form URL
#'
#' @param url Google Form URL
#' @export
#' @examples
#' \dontrun{
#' url %>% get_form_id() -> id
#' }
get_form_id <- function (url) {
    strsplit(gsub("https://docs.google.com/forms/d/", "", url),  "/")[[1]][1]
}


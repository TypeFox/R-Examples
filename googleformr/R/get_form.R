#' Scrape Google Form Content
#'
#' Get Google Form content by scraping page
#'
#' @param form Can be either the form_url or form_id
#' @export
#' @include make_url.R
#' @examples
#' \dontrun{
#' url %>% get_form() -> scrape
#' }
get_form <- function (form) {
    xml2::read_html(make_url(form, "get"))
}


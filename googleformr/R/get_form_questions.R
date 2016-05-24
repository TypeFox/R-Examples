#' Extracts Google Form Question Text
#'
#' Extracts Google Form Question Text From Scrape
#'
#' @param scrape HTML Output from Scraping Google Form
#' @export
#' @examples
#' \dontrun{
#' url %>% get_form() %>% get_form_questions() -> questions
#' }
get_form_questions <- function (scrape) {
    qt <- gsub("\n", ""
               , rvest::html_text(
                 rvest::html_nodes(scrape
                                   , ".ss-q-title") ) )
    if (length(qt) != 0) {
        qt
    } else {
      rvest::html_text(rvest::html_nodes(scrape, ".freebirdFormviewerViewItemsItemItemTitle"))
    }
}


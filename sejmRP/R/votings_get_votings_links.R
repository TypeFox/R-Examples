#' Getting votings' links
#'
#' Function \code{votings_get_votings_links} gets votings' links from
#' meeting's page.
#'
#' @details
#' Example of a meeting's page: 
#' http://www.sejm.gov.pl/Sejm7.nsf/agent.xsp?symbol=listaglos&IdDnia=1179
#'
#' @usage votings_get_votings_links(home_page = 'http://www.sejm.gov.pl/Sejm8.nsf/',
#'   page)
#'
#' @param home_page main page of polish diet: http://www.sejm.gov.pl/Sejm8.nsf/
#' @param page meeting's page
#'
#' @return character vector
#'
#' @examples
#' \dontrun{
#' home_page <- 'http://www.sejm.gov.pl/Sejm7.nsf/'
#' page <- 'http://www.sejm.gov.pl/Sejm7.nsf/agent.xsp?symbol=listaglos&IdDnia=1179'
#' votings_get_votings_links(home_page, page)}
#' 
#' @note
#' All information is stored in PostgreSQL database.
#'
#' @author Piotr Smuda
#'
#' @export
#'

votings_get_votings_links <- function(home_page = "http://www.sejm.gov.pl/Sejm8.nsf/", page) {
    stopifnot(is.character(home_page), is.character(page))
    
    # getting votings links
    votings_links <- html_nodes(safe_html(page), ".bold a")
    votings_links <- unlist(html_attrs(votings_links), use.names = FALSE)
    
    # removing votings about positions
    correct <- !stri_detect_regex(votings_links, "glosowaniaL")
    
    votings_links <- votings_links[correct]
    votings_links <- paste0(home_page, votings_links)
    
    return(votings_links)
} 

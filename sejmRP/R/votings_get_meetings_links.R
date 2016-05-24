#' Getting meetings' links
#'
#' Function \code{votings_get_meetings_links} gets meetings' links.
#'
#' @usage votings_get_meetings_links(
#'   home_page = 'http://www.sejm.gov.pl/Sejm8.nsf/', page = 
#'   'http://www.sejm.gov.pl/Sejm8.nsf/agent.xsp?symbol=posglos&NrKadencji=8')
#'
#' @param home_page main page of polish diet: http://www.sejm.gov.pl/Sejm8.nsf/
#' @param page page with votings in polish diet: 
#' http://www.sejm.gov.pl/Sejm8.nsf/agent.xsp?
#' symbol=posglos&NrKadencji=8
#'
#' @return character vector
#'
#' @examples
#' \dontrun{
#' votings_get_meetings_links()}
#' 
#' @note
#' All information is stored in PostgreSQL database.
#'
#' @author Piotr Smuda
#'
#' @export
#'

votings_get_meetings_links <- function(home_page = "http://www.sejm.gov.pl/Sejm8.nsf/",
                                       page = "http://www.sejm.gov.pl/Sejm8.nsf/agent.xsp?symbol=posglos&NrKadencji=8") {
    stopifnot(is.character(home_page), is.character(page))
    
    # getting meetings links
    meetings_links <- html_nodes(safe_html(page), "td.left")
    meetings_links <- html_nodes(meetings_links, "a")
    meetings_links <- unlist(html_attrs(meetings_links), use.names = FALSE)
    meetings_links <- paste0(home_page, meetings_links)
    
    return(meetings_links)
} 

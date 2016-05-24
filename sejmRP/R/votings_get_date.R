#' Getting date of meeting
#'
#' Function \code{votings_get_date} gets a date of meeting.
#'
#' @details
#' Example of a meeting's page: 
#' http://www.sejm.gov.pl/Sejm7.nsf/agent.xsp?symbol=listaglos&IdDnia=1179
#'
#' @usage votings_get_date(page)
#'
#' @param page meeting's page
#'
#' @return date in format YYYY-MM-DD as character
#'
#' @examples
#' \dontrun{
#' page <- 'http://www.sejm.gov.pl/Sejm7.nsf/agent.xsp?symbol=listaglos&IdDnia=1179'
#' votings_get_date(page)}
#' 
#' @note
#' All information is stored in PostgreSQL database.
#'
#' @author Piotr Smuda
#'
#' @export
#'

votings_get_date <- function(page) {
    stopifnot(is.character(page))
    
    # getting date
    date <- html_nodes(safe_html(page), "h1")
    date <- html_text(date)
    date <- unlist(stri_extract_all_regex(date, "[0-9]{2}-[0-9]{2}-[0-9]{4}"))
    date <- as.character(strptime(date, "%d-%m-%Y"))
    
    return(date)
} 

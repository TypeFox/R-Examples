#' Getting links with voting's results for each club
#'
#' Function \code{votes_get_clubs_links} gets links with voting's results for each club
#' from voting's page.
#'
#' @details
#' Function \code{votes_get_clubs_links} gets links with voting's results for each club
#' from voting's page. Example of a voting's page: 
#' http://www.sejm.gov.pl/Sejm7.nsf/agent.xsp?symbol=glosowania&
#' NrKadencji=7&NrPosiedzenia=1&NrGlosowania=1
#'
#' @usage votes_get_clubs_links(home_page = 'http://www.sejm.gov.pl/Sejm8.nsf/',
#'   page)
#'
#' @param home_page main page of polish diet: http://www.sejm.gov.pl/Sejm8.nsf/
#' @param page voting's page
#' 
#' @return data frame with two columns: club, links
#'
#' @examples
#' \dontrun{
#' home_page <- 'http://www.sejm.gov.pl/Sejm7.nsf/'
#' page <- paste0('http://www.sejm.gov.pl/Sejm7.nsf/agent.xsp?',
#'   'symbol=glosowania&NrKadencji=7&NrPosiedzenia=1&NrGlosowania=1')
#' votes_get_clubs_links(home_page, page)}
#' 
#' @note
#' All information is stored in PostgreSQL database.
#' 
#' @author Piotr Smuda
#'
#' @export
#'

votes_get_clubs_links <- function(home_page = "http://www.sejm.gov.pl/Sejm8.nsf/", page) {
    stopifnot(is.character(home_page), is.character(page))
    
    # getting clubs
    results_page <- safe_html(page)
    votes_info <- html_nodes(results_page, ".center .right")
    votes_clubs <- html_text(votes_info)
    
    # if there is comment in table
    comments <- stri_detect_regex(votes_clubs, "Uwaga|\\s")
    votes_clubs <- votes_clubs[!comments]
    
    # if there isn't table with results
    if (length(votes_clubs) == 0) {
        return(invisible(NULL))
    }
    
    # getting links
    votes_links <- html_nodes(votes_info, "a")
    votes_links <- unlist(html_attrs(votes_links), use.names = FALSE)
    votes_links <- paste0(home_page, votes_links)
    
    # creating data frame with data
    votes_clubs_links <- data.frame(club = votes_clubs, links = votes_links, stringsAsFactors = FALSE)
    
    return(votes_clubs_links)
} 

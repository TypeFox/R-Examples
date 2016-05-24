#' Getting voting's results for each club
#'
#' Function \code{votes_get_results} gets voting's results for each club.
#'
#' @details
#' Function \code{votes_get_results} gets voting's results for each club.
#' Example of page with voting's results of PO club: 
#' http://www.sejm.gov.pl/Sejm7.nsf/agent.xsp?symbol=klubglos&
#' IdGlosowania=37494&KodKlubu=PO
#'
#' @usage votes_get_results(page)
#'
#' @param page club's voting's results page
#' 
#' @return data frame with two columns: deputy, vote
#'
#' @examples
#' \dontrun{
#' page <- paste0('http://www.sejm.gov.pl/Sejm7.nsf/agent.xsp?',
#'            'symbol=klubglos&IdGlosowania=37494&KodKlubu=PO')
#' votes_get_results(page)}
#' 
#' @note
#' All information is stored in PostgreSQL database.
#' 
#' @author Piotr Smuda
#'
#' @export
#'

votes_get_results <- function(page) {
    stopifnot(is.character(page))
    
    # getting deputies and their votes
    votes_clubs_results <- safe_readHTMLTable(page, encoding = "UTF-8", stringsAsFactors = FALSE)[[1]]
    if(ncol(votes_clubs_results) == 6) {
      deputies <- c(votes_clubs_results[, 2], votes_clubs_results[, 5])
      deputies <- deputies[!is.na(deputies)]
      deputies_votes <- c(votes_clubs_results[, 3], votes_clubs_results[, 6])
      deputies_votes <- deputies_votes[!is.na(deputies_votes)]
    } else {
      deputies <- votes_clubs_results[, 2]
      deputies_votes <- votes_clubs_results[, 3]
    }
    
    # creating data frame with data
    votes_results <- data.frame(deputy = deputies, vote = deputies_votes, stringsAsFactors = FALSE)
    
    return(votes_results)
} 

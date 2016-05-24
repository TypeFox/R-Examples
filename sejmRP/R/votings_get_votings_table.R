#' Getting votings' table
#'
#' Function \code{votings_get_votings_table} gets votings' table from
#' meeting's page.
#'
#' @details
#' Function \code{votings_get_votings_table} gets votings' table from
#' meeting's page. Example of a meeting's page: 
#' http://www.sejm.gov.pl/Sejm7.nsf/agent.xsp?symbol=listaglos&IdDnia=1179
#' The result of this function is a data frame with three columns, where
#' the first includes numbers of votings, the second voting's time 
#' and the third is with voting's topics.
#'
#' @usage votings_get_votings_table(page)
#'
#' @param page meeting's page
#'
#' @return data frame with three columns: Nr, Godzina (Time), Temat (Topic)
#'
#' @examples
#' \dontrun{
#' page <- 'http://www.sejm.gov.pl/Sejm7.nsf/agent.xsp?symbol=listaglos&IdDnia=1179'
#' votings_get_votings_table(page)}
#' 
#' @note
#' All information is stored in PostgreSQL database.
#'
#' @author Piotr Smuda
#'
#' @export
#'

votings_get_votings_table <- function(page) {
    stopifnot(is.character(page))
    
    # getting votings table
    votings_table <- safe_readHTMLTable(page, encoding = "UTF-8", stringsAsFactors = FALSE)[[1]]
    
    return(votings_table)
} 

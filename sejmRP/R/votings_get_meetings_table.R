#' Getting meetings' table
#'
#' Function \code{votings_get_meetings_table} gets meetings' table.
#'
#' @details
#' Function \code{votings_get_meetings_table} gets meetings' table. The
#' result of this function is a data frame with three columns, where
#' the first includes numbers of meetings, the second theirs' dates in
#' Polish and the third is with numbers of votings on each meeting.
#'
#' @usage votings_get_meetings_table(page = 
#'   'http://www.sejm.gov.pl/Sejm8.nsf/agent.xsp?symbol=posglos&NrKadencji=8')
#'
#' @param page page with votings in polish diet: 
#' http://www.sejm.gov.pl/Sejm8.nsf/agent.xsp?
#' symbol=posglos&NrKadencji=8
#'
#' @return data frame with three unnamed columns
#'
#' @examples
#' \dontrun{
#' votings_get_meetings_table()}
#' 
#' @note
#' All information is stored in PostgreSQL database.
#'
#' @author Piotr Smuda
#'
#' @export
#'

votings_get_meetings_table <- function(page = "http://www.sejm.gov.pl/Sejm8.nsf/agent.xsp?symbol=posglos&NrKadencji=8") {
    stopifnot(is.character(page))
    
    # getting meetings table
    meetings_table <- safe_readHTMLTable(page, encoding = "UTF-8", stringsAsFactors = FALSE)[[1]]
    
    # filling first column where number of meeting is missing
    meeting_number <- meetings_table[1, 1]
    for (i in seq_len(nrow(meetings_table))) {
        if (meetings_table[i, 1] != "") {
            meeting_number <- meetings_table[i, 1]
        }
        meetings_table[i, 1] <- meeting_number
    }
    
    return(meetings_table)
} 

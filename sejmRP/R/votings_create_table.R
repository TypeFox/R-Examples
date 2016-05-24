#' Creating table with votings
#'
#' Function \code{votings_create_table} creates a table with votings.
#' 
#' @usage votings_create_table(dbname, user, password, host,
#'   nr_term_of_office = 8)
#'
#' @param dbname name of database
#' @param user name of user
#' @param password password of database
#' @param host name of host
#' @param nr_term_of_office number of term of office of Polish Diet; default: 8
#'
#' @return invisible NULL
#'
#' @examples
#' \dontrun{
#' votings_create_table(dbname, user, password, host)}
#' 
#' @note
#' Use only this function for first time, when the \emph{votings} table
#' is empty. Then use \code{votings_update_table}.
#' 
#' All information is stored in PostgreSQL database.
#'
#' @author Piotr Smuda
#'
#' @export
#'

votings_create_table <- function(dbname, user, password, host, nr_term_of_office = 8) {
    stopifnot(is.character(dbname), is.character(user), is.character(password), is.character(host), 
              is.numeric(nr_term_of_office), nr_term_of_office%%1 == 0)
    
    # setting home page and page with votings
    home_page <- paste0("http://www.sejm.gov.pl/Sejm", nr_term_of_office, ".nsf/")
    page <- paste0("http://www.sejm.gov.pl/Sejm", nr_term_of_office,
                   ".nsf/agent.xsp?symbol=posglos&NrKadencji=", nr_term_of_office)
  
    # getting meetings table with meetings' numbers
    meetings_table <- votings_get_meetings_table(page)
    
    # getting meetings links with votings
    meetings_links <- votings_get_meetings_links(home_page, page)
    
    # remembering first voting id
    id_voting <- 1
    
    for (i in rev(seq_len(length(meetings_links)))) {
        # getting meetings date
        meetings_date <- votings_get_date(meetings_links[i])
        
        # getting votings table with votings' numbers and topics
        votings_table <- votings_get_votings_table(meetings_links[i])
        
        # getting votings links
        votings_links <- votings_get_votings_links(home_page, meetings_links[i])
        
        # putting this data frame to database
        drv <- dbDriver("PostgreSQL")
        database_diet <- dbConnect(drv, dbname = dbname, user = user, password = password, host = host)
        for (j in rev(seq_len(length(votings_links)))) {
            dbSendQuery(database_diet, paste0("INSERT INTO votings (id_voting, nr_term_of_office, nr_meeting, date_meeting, ",
                                              "nr_voting, topic_voting, link_results) VALUES (", id_voting, ",",
                                              nr_term_of_office, ",", meetings_table[i, 1], ",'", meetings_date, "',", 
                                              votings_table[j, 1], ",'", votings_table[j, 3], "','", votings_links[j], "')"))
            
            id_voting <- id_voting + 1
        }
        suppressWarnings(dbDisconnect(database_diet))
    }
    
    # creating a flag file when creating votings table is finished
    file.create("sejmRP_votings_flag")
    
    return(invisible(NULL))
} 

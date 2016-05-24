#' Matching deputies to theirs' ids
#'
#' Function \code{votes_match_deputies_ids} matches deputies from voting's results
#' page to theirs' ids from \emph{deputies} table.
#'
#' @details
#' Function \code{votes_match_deputies_ids} matches deputies from voting's results
#' page to theirs' ids from \emph{deputies} table. The result of this function is
#' a data frame with deputies' data, ids and votes. Because of encoding issue
#' on Windows operation system, you need to select if you use Windows.
#' Example of page with voting's results of PO club: 
#' http://www.sejm.gov.pl/Sejm7.nsf/agent.xsp?
#' symbol=klubglos&IdGlosowania=37494&KodKlubu=PO
#'
#' @usage votes_match_deputies_ids(dbname, user, password, host, page, 
#'    nr_term_of_office = 8, windows = .Platform$OS.type == 'windows')
#'
#' @param dbname name of database
#' @param user name of user
#' @param password password of database
#' @param host name of host
#' @param page club's voting's results page
#' @param nr_term_of_office number of term of office of Polish Diet; default: 8
#' @param windows information of used operation system; 
#' default: .Platform$OS.type == 'windows'
#' 
#' @return data frame with three columns: deputy, vote, id
#'
#' @examples
#' \dontrun{
#' page <- paste0('http://www.sejm.gov.pl/Sejm7.nsf/agent.xsp?',
#'              'symbol=klubglos&IdGlosowania=37494&KodKlubu=PO')
#' votes_match_deputies_ids(dbname, user, password, host, page, 7, TRUE)
#' votes_match_deputies_ids(dbname, user, password, host, page, 7, FALSE)}
#' 
#' @note
#' All information is stored in PostgreSQL database.
#' 
#' @author Piotr Smuda
#'
#' @export
#'

votes_match_deputies_ids <- function(dbname, user, password, host, page, nr_term_of_office = 8,
                                     windows = .Platform$OS.type == "windows") {
    stopifnot(is.character(dbname), is.character(user), is.character(password), is.character(host), 
              is.character(page), is.numeric(nr_term_of_office), nr_term_of_office%%1 == 0,
              is.logical(windows))
    
    # getting all of deputies' ids
    deputies_whole_ids <- deputies_get_ids(dbname, user, password, host, nr_term_of_office, windows)
    
    # getting votes' results
    votes_results <- votes_get_results(page)
    
    # getting deputies' (from voting) ids
    deputies_ids <- unname(deputies_whole_ids[votes_results[, 1]])
    
    # creating data frame with data
    votes_deputies_ids <- cbind(votes_results, id = deputies_ids)
    
    return(votes_deputies_ids)
} 

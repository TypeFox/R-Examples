#' Updating table with votes
#'
#' Function \code{votes_update_table} updates a table with votes.
#'
#' @usage votes_update_table(dbname, user, password, host,
#'   nr_term_of_office = 8, windows = .Platform$OS.type == 'windows',
#'   verbose = FALSE)
#'
#' @param dbname name of database
#' @param user name of user
#' @param password password of database
#' @param host name of host
#' @param nr_term_of_office number of term of office of Polish Diet; default: 8
#' @param windows information of used operation system; 
#' default: .Platform$OS.type == 'windows'
#' @param verbose if TRUE then additional info will be printed
#' 
#' @return invisible NULL
#'
#' @examples
#' \dontrun{
#' votes_update_table(dbname, user, password, host, 7, TRUE)
#' votes_update_table(dbname, user, password, host, 7, FALSE)}
#' 
#' @note
#' There is a possibility that someone's voice reader broke during voting
#' and this situation is treated like this deputy was absent. Even if deputy
#' made a decision, he's/she's vote is 'Nieobecny'.
#' 
#' All information is stored in PostgreSQL database.
#' 
#' @author Piotr Smuda
#'
#' @export
#'

votes_update_table <- function(dbname, user, password, host, nr_term_of_office = 8, 
                               windows = .Platform$OS.type == "windows",
                               verbose = FALSE) {
    stopifnot(is.character(dbname), is.character(user), is.character(password), is.character(host),
              is.numeric(nr_term_of_office), nr_term_of_office%%1 == 0, is.logical(windows),
              is.logical(verbose))
    
    # setting home page and page with votings
    home_page <- paste0("http://www.sejm.gov.pl/Sejm", nr_term_of_office, ".nsf/")
  
    # checking last nr_meeting, removing records with that number and checking last id_voting
    drv <- dbDriver("PostgreSQL")
    database_diet <- dbConnect(drv, dbname = dbname, user = user, password = password, host = host)
    last_id_voting <- dbGetQuery(database_diet, paste0("SELECT max(id_voting) FROM votes WHERE nr_term_of_office = ",
                                                       nr_term_of_office))
    last_id_voting <- as.integer(last_id_voting)
    if (is.na(last_id_voting)) {
      last_id_voting <- 1
    }
    dbSendQuery(database_diet, paste0("DELETE FROM votes WHERE nr_term_of_office = ",
                                      nr_term_of_office, " and id_voting >= ", last_id_voting))
    last_id_vote <- dbGetQuery(database_diet, "SELECT max(id_vote) FROM votes")
    last_id_vote <- as.integer(last_id_vote)
    
    # getting voting_id and results links
    votings_ids_links <- dbGetQuery(database_diet, paste0("SELECT * FROM votings WHERE nr_term_of_office = ",
                                                          nr_term_of_office, " ORDER BY id_voting"))[, c(1, 7)]
    suppressWarnings(dbDisconnect(database_diet))
    
    # choosing the place from where we start update table
    which_to_update <- which(votings_ids_links[, 1] >= last_id_voting)
    votings_ids_links <- votings_ids_links[which_to_update, ]
    
    # remembering first new vote id
    id_vote <- last_id_vote + 1
    
    for (i in seq_len(nrow(votings_ids_links))) {
        # getting clubs and links from voting
      if (verbose) {
        cat("Downloading",votings_ids_links[i, 2],"\n")
      }
      votes_get_clubs <- votes_get_clubs_links(home_page, votings_ids_links[i, 2])
        
        # if there isn't table with results
        if (is.null(votes_get_clubs)) {
            next
        }
        
        for (j in seq_len(nrow(votes_get_clubs))) {
            # getting deputies id, vote and club
            votes_info <- votes_match_deputies_ids(dbname, user, password, host, votes_get_clubs[j, 2],
                                                   nr_term_of_office, windows)
            
            # putting this data frame to database
            database_diet <- dbConnect(drv, dbname = dbname, user = user, password = password, host = host)
            for (k in seq_len(nrow(votes_info))) {
                dbSendQuery(database_diet, paste0("INSERT INTO votes (id_vote, nr_term_of_office, id_deputy,",
                                                  "id_voting, vote, club) VALUES (", id_vote, ",",
                                                  nr_term_of_office, ",'", votes_info[k, 3], "',", 
                                                  votings_ids_links[i, 1], ",'", votes_info[k, 2], "','", 
                                                  votes_get_clubs[j, 1], "')"))
                id_vote <- id_vote + 1
            }
            suppressWarnings(dbDisconnect(database_diet))
        }
    }
    
    return(invisible(NULL))
} 

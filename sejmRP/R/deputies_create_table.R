#' Creating table with deputies
#'
#' Function \code{deputies_create_table} creates a table with deputies.
#'
#' @usage deputies_create_table(dbname, user, password, host,
#'    nr_term_of_office = 8)
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
#' deputies_create_table(dbname, user, password, host)}
#'
#' @note
#' Use only this function for first time, when the \emph{deputies} table
#' is empty. Then use \code{deputies_update_table}.
#' 
#' All information is stored in PostgreSQL database.
#'
#' @author Piotr Smuda
#'
#' @export
#'

deputies_create_table <- function(dbname, user, password, host, nr_term_of_office = 8) {
    stopifnot(is.character(dbname), is.character(user), is.character(password), is.character(host),
              is.numeric(nr_term_of_office), nr_term_of_office%%1 == 0)
    
    # getting data from page with active deputies
    deputies_active <- deputies_get_data("active", nr_term_of_office)
    
    # getting data from page with inactive deputies
    deputies_inactive <- deputies_get_data("inactive", nr_term_of_office)
    
    # merging ids, names and surnames of active and inactive deputies
    id_deputies <- c(deputies_active[, 1], deputies_inactive[, 1])
    deputies <- c(deputies_active[, 2], deputies_inactive[, 2])
    
    # putting this data frame to database
    drv <- dbDriver("PostgreSQL")
    database_diet <- dbConnect(drv, dbname = dbname, user = user, password = password, host = host)
    n <- length(deputies)
    for (i in seq_len(n)) {
        dbSendQuery(database_diet, paste0("INSERT INTO deputies (id_deputy, nr_term_of_office, surname_name) VALUES ('",
            id_deputies[i], "',", nr_term_of_office, ",'", deputies[i], "')"))
    }
    suppressWarnings(dbDisconnect(database_diet))
    return(invisible(NULL))
} 

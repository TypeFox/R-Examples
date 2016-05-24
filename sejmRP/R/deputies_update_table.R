#' Updating table with deputies
#'
#' Function \code{deputies_update_table} updates a table with deputies.
#'
#' @usage deputies_update_table(dbname, user, password, host,
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
#' deputies_update_table(dbname, user, password, host)}
#'
#' @note
#' All information is stored in PostgreSQL database.
#'
#' @author Piotr Smuda
#'
#' @export
#'

deputies_update_table <- function(dbname, user, password, host, nr_term_of_office = 8) {
    stopifnot(is.character(dbname), is.character(user), is.character(password), is.character(host),
              is.numeric(nr_term_of_office), nr_term_of_office%%1 == 0)
    
    # checking last id of deputies
    drv <- dbDriver("PostgreSQL")
    database_diet <- dbConnect(drv, dbname = dbname, user = user, password = password, host = host)
    last_id <- dbGetQuery(database_diet, paste0("SELECT max(id_deputy) FROM deputies WHERE nr_term_of_office = ",
                                                nr_term_of_office))
    last_id <- as.character(last_id)
    if (is.na(last_id)) {
        last_id <- "000"
    }
    suppressWarnings(dbDisconnect(database_diet))
    
    # adding new deputies to database
    deputies_add_new(dbname, user, password, host, "active", last_id, nr_term_of_office)
    deputies_add_new(dbname, user, password, host, "inactive", last_id, nr_term_of_office)
    
    return(invisible(NULL))
} 

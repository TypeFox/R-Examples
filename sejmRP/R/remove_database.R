#' Removing database
#'
#' Function \code{remove_database} remove whole database.
#'
#' @usage remove_database(dbname, user, password, host)
#'
#' @param dbname name of database
#' @param user name of user
#' @param password password of database
#' @param host name of host
#'
#' @return invisible NULL
#'
#' @examples
#' \dontrun{
#' remove_database(dbname, user, password, host)}
#' 
#' @note
#' All information is stored in PostgreSQL database.
#'
#' @author Piotr Smuda
#'
#' @export
#'

remove_database <- function(dbname, user, password, host) {
    stopifnot(is.character(dbname), is.character(user), is.character(password), is.character(host))
    
    # connecting to database
    drv <- dbDriver("PostgreSQL")
    database_diet <- dbConnect(drv, dbname = dbname, user = user, password = password, host = host)
    
    # removing tables
    dbRemoveTable(database_diet, "statements")
    dbRemoveTable(database_diet, "votes")
    dbRemoveTable(database_diet, "votings")
    dbRemoveTable(database_diet, "deputies")
    
    # disconnecting to database
    suppressWarnings(dbDisconnect(database_diet))
    return(invisible(NULL))
} 

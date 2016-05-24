#' Getting deputies' ids
#'
#' Function \code{deputies_get_ids} gets deputies' ids from \emph{deputies} table.
#'
#' @details
#' Function \code{deputies_get_ids} gets deputies' ids from \emph{deputies} table.
#' As result of this function you get named character vector with ids, where their
#' names are names and surnames of deputies. Because of encoding issue on Windows
#' operation system, you need to select if you use Windows.
#'
#' @usage deputies_get_ids(dbname, user, password, host,
#'    nr_term_of_office = 8, windows = .Platform$OS.type == 'windows')
#'
#' @param dbname name of database
#' @param user name of user
#' @param password password of database
#' @param host name of host
#' @param nr_term_of_office number of term of office of Polish Diet; default: 8
#' @param windows information of used operation system; default: 
#' .Platform$OS.type == 'windows'
#' 
#' @return named character vector
#'
#' @examples
#' \dontrun{
#' deputies_get_ids(dbname, user, password, host, TRUE)
#' deputies_get_ids(dbname, user, password, host, FALSE)}
#' 
#' @note
#' All information is stored in PostgreSQL database.
#'
#' @author Piotr Smuda
#'
#' @export
#'

deputies_get_ids <- function(dbname, user, password, host, nr_term_of_office = 8, windows = .Platform$OS.type == "windows") {
    stopifnot(is.character(dbname), is.character(user), is.character(password), is.character(host),
              is.numeric(nr_term_of_office), nr_term_of_office%%1 == 0, is.logical(windows))
    
    # getting deputies' ids
    drv <- dbDriver("PostgreSQL")
    database_diet <- dbConnect(drv, dbname = dbname, user = user, password = password, host = host)
    deputies_table <- dbGetQuery(database_diet, paste0("SELECT * FROM deputies WHERE nr_term_of_office = ", nr_term_of_office))
    deputies_ids <- deputies_table[, 1]
    
    # calling deputies' ids with their names and surnames
    if (windows) {
        names(deputies_ids) <- iconv(deputies_table[, 3], from = "UTF-8", to = "Windows-1250")
    } else {
        names(deputies_ids) <- deputies_table[, 3]
    }
    suppressWarnings(dbDisconnect(database_diet))
    
    return(deputies_ids)
} 

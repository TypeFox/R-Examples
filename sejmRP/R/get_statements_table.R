#' Importing statements table from a database
#'
#' Function \code{get_statements_table} imports statements table from a database.
#'
#' @details
#' Function \code{get_statements_table} imports statements table from a database.
#' The result of this function is a data frame with statements' data. Because of
#' encoding issue on Windows operation system, you need to select if you use Windows.
#'
#' @usage get_statements_table(dbname = 'sejmrp', user = 'reader',
#'   password = 'qux94874', host = 'services.mini.pw.edu.pl',
#'   sorted_by_id = TRUE, windows = .Platform$OS.type == 'windows')
#'
#' @param dbname name of database; default: 'sejmrp'
#' @param user name of user; default: 'reader'
#' @param password password of database; default: 'qux94874'
#' @param host name of host; default: 'services.mini.pw.edu.pl'
#' @param sorted_by_id information if table should be sorted by id; default: TRUE
#' @param windows information of used operation system;
#' default: .Platform$OS.type == 'windows'
#'
#' @return data frame
#'
#' @examples
#' \dontrun{
#' statements <- get_statements_table()
#' dim(statements)
#' # [1] 43432 6
#' names(statements)
#' # [1] 'id_statement' 'nr_term_of_office' 'surname_name'
#' # [4] 'date_statement' 'titles_order_points' 'statement'}
#'
#' @note
#' Default parameters use privilages of 'reader'. It can only SELECT data from database.
#'
#' All information is stored in PostgreSQL database.
#'
#' @author Piotr Smuda
#'
#' @export
#'

get_statements_table <- function(dbname = "sejmrp", user = "reader", password = "qux94874", host = "services.mini.pw.edu.pl",
    sorted_by_id = TRUE, windows = .Platform$OS.type == "windows") {
    stopifnot(is.character(dbname), is.character(user), is.character(password), is.character(host), is.logical(sorted_by_id),
        is.logical(windows))

    # connecting to database
    drv <- dbDriver("PostgreSQL")
    database_diet <- dbConnect(drv, dbname = dbname, user = user, password = password, host = host)

    # add information about new SELECT to the counter table
    dbSendQuery(database_diet, paste0("INSERT INTO counter (what, date) VALUES ('statements','", Sys.Date(), "')"))

    # reading table
    if (sorted_by_id) {
        statements <- dbGetQuery(database_diet, "SELECT * FROM statements ORDER BY nr_term_of_office, id_statement")
    } else {
        statements <- dbGetQuery(database_diet, "SELECT * FROM statements")
    }

    # encoding for windows
    if (windows) {
        statements[, 3] <- iconv(statements[, 3], from = "UTF-8", to = "Windows-1250")
        statements[, 5] <- iconv(statements[, 5], from = "UTF-8", to = "Windows-1250")
        statements[, 6] <- iconv(statements[, 6], from = "UTF-8", to = "Windows-1250")
    }

    suppressWarnings(dbDisconnect(database_diet))
    return(statements)
}

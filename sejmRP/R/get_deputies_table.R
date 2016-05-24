#' Importing deputies table from a database
#'
#' Function \code{get_deputies_table} imports deputies table from a database.
#'
#' @details
#' Function \code{get_deputies_table} imports deputies table from a database.
#' The result of this function is a data frame with deputies' data. Because of
#' encoding issue on Windows operation system, you need to select if you use Windows.
#'
#' @usage get_deputies_table(dbname = 'sejmrp', user = 'reader',
#'   password = 'qux94874', host = 'services.mini.pw.edu.pl',
#'   sorted_by_id = TRUE, windows = .Platform$OS.type == 'windows')
#'
#' @param dbname name of database; default: 'sejmrp'
#' @param user name of user; default: 'reader'
#' @param password password of database; default: 'qux94874'
#' @param host name of host; default: 'services.mini.pw.edu.pl'
#' @param sorted_by_id information if table should be sorted by id; default: TRUE
#' @param windows information of used operation system; default: .Platform$OS.type == 'windows'
#'
#' @return data frame
#'
#' @examples
#' \dontrun{
#' deputies <- get_deputies_table()
#' dim(deputies)
#' # [1] 983 3
#' names(deputies)
#' # [1] 'id_deputy' 'nr_term_of_office' 'surname_name'}
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

get_deputies_table <- function(dbname = "sejmrp", user = "reader", password = "qux94874", host = "services.mini.pw.edu.pl",
                               sorted_by_id = TRUE, windows = .Platform$OS.type == "windows") {
    stopifnot(is.character(dbname), is.character(user), is.character(password), is.character(host), is.logical(sorted_by_id),
        is.logical(windows))

    # connecting to database
    drv <- dbDriver("PostgreSQL")
    database_diet <- dbConnect(drv, dbname = dbname, user = user, password = password, host = host)

    # add information about new SELECT to the counter table
    dbSendQuery(database_diet, paste0("INSERT INTO counter (what, date) VALUES ('deputies','", Sys.Date(), "')"))

    # reading table
    if (sorted_by_id) {
        deputies <- dbGetQuery(database_diet, "SELECT * FROM deputies ORDER BY nr_term_of_office, id_deputy")
    } else {
        deputies <- dbGetQuery(database_diet, "SELECT * FROM deputies")
    }

    # encoding for windows
    if (windows) {
        deputies[, 3] <- iconv(deputies[, 3], from = "UTF-8", to = "Windows-1250")
    }

    suppressWarnings(dbDisconnect(database_diet))
    return(deputies)
}

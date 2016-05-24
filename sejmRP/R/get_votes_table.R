#' Importing votes table from a database
#'
#' Function \code{get_votes_table} imports votes table from a database.
#'
#' @details
#' Function \code{get_votes_table} imports votes table from a database.
#' The result of this function is a data frame with votes' data. Because of
#' encoding issue on Windows operation system, you need to select if you use Windows.
#'
#' @usage get_votes_table(dbname = 'sejmrp', user = 'reader',
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
#' votes <- get_votes_table()
#' dim(votes)
#' # [1] 2826483 6
#' names(votes)
#' # [1] 'id_vote' 'nr_term_of_office' 'id_deputy' 'id_voting' 'vote' 'club'
#' object.size(votes)
#' # 90474040 bytes}
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

get_votes_table <- function(dbname = "sejmrp", user = "reader", password = "qux94874", host = "services.mini.pw.edu.pl",
                            sorted_by_id = TRUE, windows = .Platform$OS.type == "windows") {
    stopifnot(is.character(dbname), is.character(user), is.character(password), is.character(host),
              is.logical(sorted_by_id), is.logical(windows))

    # connecting to database
    drv <- dbDriver("PostgreSQL")
    database_diet <- dbConnect(drv, dbname = dbname, user = user, password = password, host = host)

    # add information about new SELECT to the counter table
    dbSendQuery(database_diet, paste0("INSERT INTO counter (what, date) VALUES ('votes','", Sys.Date(), "')"))

    # reading table
    if (sorted_by_id) {
        votes <- dbGetQuery(database_diet, "SELECT * FROM votes ORDER BY nr_term_of_office, id_voting, id_vote")
    } else {
        votes <- dbGetQuery(database_diet, "SELECT * FROM votes")
    }

    # encoding for windows
    if (windows) {
        votes[, 5] <- iconv(votes[, 5], from = "UTF-8", to = "Windows-1250")
        votes[, 6] <- iconv(votes[, 6], from = "UTF-8", to = "Windows-1250")
    }

    suppressWarnings(dbDisconnect(database_diet))
    return(votes)
}

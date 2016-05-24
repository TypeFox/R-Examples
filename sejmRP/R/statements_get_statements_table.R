#' Getting statements' table
#'
#' Function \code{statements_get_statements_table} gets statements' table from
#' meeting's page.
#'
#' @details
#' Function \code{statements_get_statements_table} gets statements' table. from
#' meeting's page. Example of a meeting's page: 
#' http://www.sejm.gov.pl/Sejm7.nsf/posiedzenie.xsp?posiedzenie=99&dzien=2
#' The result of this function is a data frame with three columns, where
#' the first includes author of statement, the second the number of order point 
#' and the third is a title of order point.
#'
#' @usage statements_get_statements_table(page)
#'
#' @param page meeting's page
#'
#' @return data frame with three unnamed columns
#'
#' @examples
#' \dontrun{
#' page <- 'http://www.sejm.gov.pl/Sejm7.nsf/posiedzenie.xsp?posiedzenie=99&dzien=2'
#' statements_get_statements_table(page)}
#' 
#' @note
#' All information is stored in PostgreSQL database.
#'
#' @author Piotr Smuda
#'
#' @export
#'

statements_get_statements_table <- function(page) {
    stopifnot(is.character(page))
        
    # getting statements' table
    statements_table <- safe_readHTMLTable(page, encoding = "UTF-8", stringsAsFactors = FALSE)[[1]]
    statements_table[, 3] <- stri_replace_all_regex(statements_table[, 3], "'", "")
    
    # changing colnames to numbers
    names(statements_table) <- seq_len(ncol(statements_table))
    
    return(statements_table)
} 
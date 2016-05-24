#' Getting statements
#'
#' Function \code{statements_get_statement} gets statement's content.
#' 
#' @details
#' Function \code{statements_get_statement} gets statement's content.
#' Example of page with deputy's statement: 
#' http://www.sejm.gov.pl/Sejm7.nsf/wypowiedz.xsp?posiedzenie=15&dzien=1&wyp=008
#'
#' @usage statements_get_statement(page, ...)
#'
#' @param page deputy's statement's page
#' @param ... other arguments, that will be passed to safe_html()
#'
#' @return character vector
#'
#' @examples
#' \dontrun{
#' page <- paste0('http://www.sejm.gov.pl/Sejm7.nsf/',
#'                'wypowiedz.xsp?posiedzenie=15&dzien=1&wyp=008')
#' statements_get_statement(page)}
#'
#' @note
#' All information is stored in PostgreSQL database.
#'
#' @author Piotr Smuda, Tomasz Mikolajczyk
#'
#' @export
#'

statements_get_statement <- function(page, ...) {
    stopifnot(is.character(page))
    
    pageH <- safe_html(page, ...)
    page <- html_nodes(pageH, ".stenogram p")
    
    # getting statement content
    statement <- html_text(page)
    statement <- stri_trim_both(statement)
    statement <- paste0(statement, collapse = " ")
    statement <- stri_replace_all_regex(statement, "'", "")
    
    return(statement)
} 

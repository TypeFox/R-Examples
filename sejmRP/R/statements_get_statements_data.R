#' Getting data about statements
#'
#' Function \code{statements_get_statements_data} gets data about statements.
#' 
#' @details
#' Function \code{statements_get_statements_data} gets data about statements like
#' author, page with content of statement and it's id.
#'
#' @usage statements_get_statements_data(statements_links,
#'   home_page = 'http://www.sejm.gov.pl/')
#'
#' @param statements_links list of elements of XMLNodeSet class with statements' ids, links 
#' and their's authors
#' @param home_page main page of polish diet: http://www.sejm.gov.pl/
#'
#' @return data frame with three columns: names, statements_links, ids
#'
#' @examples
#' \dontrun{
#' page <- safe_html(paste0('http://www.sejm.gov.pl/Sejm7.nsf/',
#'                'wypowiedz.xsp?posiedzenie=15&dzien=1&wyp=0'))
#' page <- html_nodes(page, '.stenogram')
#' statements_links <- html_nodes(page, 'h2 a')
#' statements_get_statements_data(statements_links,
#'   home_page = 'http://www.sejm.gov.pl/Sejm7.nsf/')}
#'
#' @note
#' All information is stored in PostgreSQL database.
#'
#' @author Piotr Smuda, Tomasz Mikolajczyk
#'
#' @export
#'

statements_get_statements_data <- function(statements_links, home_page = "http://www.sejm.gov.pl/") {
    stopifnot(is.list(statements_links), length(statements_links) > 0)
    
    # Get speakers' names
    names <- html_text(statements_links)
    
    # Check whether the speaker is a deputy
    if_deputy <- stri_detect_regex(names, "(Pose.{1,2} )|(Minister )|([p|P]rezes Rady Ministr.{1,2} )")
    names <- names[if_deputy]
    
    # Get links to statements
    statements_links <- html_attr(statements_links, "href")
    statements_links <- statements_links[if_deputy]
    statements_links <- paste0(home_page, statements_links)
    
    # Get statements' ids
    ids <- unlist(stri_extract_all_regex(statements_links, "[0-9]+$"))
    
    # Create dataframe with deputies' names and statements and statements' ids
    statements_data <- data.frame(names = names, statements_links = statements_links, ids = ids, stringsAsFactors = FALSE)
    
    return(statements_data)
} 

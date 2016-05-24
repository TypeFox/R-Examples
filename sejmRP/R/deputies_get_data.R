#' Getting data about deputies
#'
#' Function \code{deputies_get_data} gets data about deputies.
#'
#' @details
#' Function \code{deputies_get_data} gets deputies' ids and personal data like
#' name and surname. Also there is a choice between types of deputies, because
#' on the page of Polish diet deputies are splitted into \emph{active} and \emph{inactive}.
#'
#' @usage deputies_get_data(type, nr_term_of_office = 8)
#'
#' @param type type of deputies which be add to table with deputies: active, inactive
#' @param nr_term_of_office number of term of office of Polish Diet; default: 8
#'
#' @return data frame with two columns: id_deputy, surname_name
#'
#' @examples
#' \dontrun{
#' deputies_get_data('active')
#' deputies_get_data('inactive')}
#'
#' @note
#' All information is stored in PostgreSQL database.
#'
#' @author Piotr Smuda
#'
#' @export
#'

deputies_get_data <- function(type, nr_term_of_office = 8) {
    
    stopifnot(is.character(type), type == "active" || type == "inactive",
              is.numeric(nr_term_of_office), nr_term_of_office%%1 == 0)
    
    # choosing proper page of deputies
    if (type == "active") {
        page <- paste0("http://www.sejm.gov.pl/Sejm", nr_term_of_office, ".nsf/poslowie.xsp?type=A")
    } else if (type == "inactive") {
        page <- paste0("http://www.sejm.gov.pl/Sejm", nr_term_of_office, ".nsf/poslowie.xsp?type=B")
    }
    
    # getting data from page with deputies
    page <- safe_html(page)
    deputies_data <- html_nodes(page, "#contentBody li")
    
    # getting ids of deputies
    links <- html_nodes(deputies_data, "a")
    links <- html_attrs(links)
    links <- sapply(links, function(element) {
        element[2]
    })
    id_deputies <- unlist(stri_extract_all_regex(links, "(?<=id=)[0-9]+"))
    
    # getting names and surnnames of deputies
    deputies <- html_nodes(deputies_data, ".deputyName")
    deputies <- html_text(deputies)
    
    # creating data frame with deputies data
    deputies_df <- data.frame(id_deputy = id_deputies, surname_name = deputies, stringsAsFactors = FALSE)
    
    return(deputies_df)
} 

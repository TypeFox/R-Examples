#' Creating table with deputies' statements
#'
#' Function \code{statements_create_table} creates a table with deputies' statements.
#'
#' @usage statements_create_table(dbname, user, password, host,
#'   nr_term_of_office = 8)
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
#' statements_create_table(dbname, user, password, host)}
#'
#' @note
#' Use only this function for first time, when the \emph{statements} table
#' is empty. Then use \code{statements_update_table}.
#' 
#' All information is stored in PostgreSQL database.
#'
#' @author Piotr Smuda, Tomasz Mikolajczyk
#'
#' @export
#'

statements_create_table <- function(dbname, user, password, host, nr_term_of_office = 8) {
    stopifnot(is.character(dbname), is.character(user), is.character(password), is.character(host), 
              is.numeric(nr_term_of_office), nr_term_of_office%%1 == 0)
    
    # set home page and pattern for css in statements_links
    if (nr_term_of_office == 7) {
      home_page <- home_page_links <- "http://www.sejm.gov.pl/Sejm7.nsf/"
      css_pattern <- "h2 a"
    } else if (nr_term_of_office == 8){
      home_page_links <- "http://www.sejm.gov.pl"
      home_page <- paste0(home_page_links, "/Sejm8.nsf/")
      css_pattern <- ".mowca-link a"
    }
  
    # set meeting and day number to 1
    nr_meeting <- 1
    nr_day <- 1
    repeat {
        repeat {
            # get statements links of first day of a meeting
            page <- paste0(home_page, "wypowiedz.xsp?posiedzenie=", nr_meeting, "&dzien=", nr_day, "&wyp=0")
            stenogram <- html_nodes(safe_html(page), ".stenogram")
            statements_links <- html_nodes(stenogram, css_pattern)
            
            # move to next day of meeting if empty page found
            if (length(statements_links) == 0) {
                break
            }
          
            # get titles of order points during a meeting
            page_meeting <- paste0(home_page, "posiedzenie.xsp?posiedzenie=", nr_meeting, "&dzien=", nr_day)
            statements_table <- statements_get_statements_table(page_meeting)
            if_deputy <- stri_detect_regex(statements_table[, 1], "(Pose.{1,2} )|(Minister )|([p|P]rezes Rady Ministr.{1,2} )")
            titles_order_points <- statements_table[if_deputy, 3]
            
            # get date
            statements_date <- votings_get_date(page)
            
            # get deputies' names, statements and statements' ids
            statements_data <- statements_get_statements_data(statements_links, home_page_links)
            
            # get statements
            statements <- unlist(lapply(statements_data[, 2], function(elem) {
                statement <- statements_get_statement(elem)
            }))
            
            # put dataframe to database
            drv <- dbDriver("PostgreSQL")
            database_diet <- dbConnect(drv, dbname = dbname, user = user, password = password, host = host)
            
            for (i in seq_len(length(statements))) {
                id <- paste0(nr_meeting, ".", nr_day, ".", statements_data[i, 3])
                dbSendQuery(database_diet, paste0("INSERT INTO statements (id_statement, nr_term_of_office, surname_name, date_statement, ", 
                  "titles_order_points, statement) VALUES ('", id, "',", nr_term_of_office, ",'", statements_data[i, 1], "','",
                  statements_date, "','", titles_order_points[i], "','", statements[i], "')"))
            }
            
            suppressWarnings(dbDisconnect(database_diet))
            
            # next day of meeting
            nr_day <- nr_day + 1
        }
        
        # break if last meeting found
        if (nr_day == 1) {
            break
        }
        # next meeting
        nr_day <- 1
        nr_meeting <- nr_meeting + 1
    }
    
    return(invisible(NULL))
} 

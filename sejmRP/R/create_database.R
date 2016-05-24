#' Creating database
#'
#' Function \code{create_database} creates a database with four empty
#' tables: deputies, votings, votes, statements.
#' 
#' @details
#' \preformatted{
#' Created tables:
#' 1. deputies with columns:
#'     1) id_deputy - deputy's id,
#'     2) nr_term_of_office - Polish Diet's number of term of office,
#'     3) surname_name - deputy's names and surnames,
#' 2. votings with columns:
#'     1) id_voting - voting's id,
#'     2) nr_term_of_office - Polish Diet's number of term of office,
#'     3) nr_meeting - meeting's number,
#'     4) date_meeting - meeting's date,
#'     5) nr_voting - voting's number,
#'     6) topic_voting - voting's topic,
#'     7) link_results - link with voting's results,
#' 3. votes with columns:
#'     1) id_vote - vote's id,
#'     2) nr_term_of_office - Polish Diet's number of term of office,
#'     3) id_deputy - deputy's id,
#'     4) id_voting - voting's id,
#'     5) vote - deputy's vote, one of: 'Za','Przeciw',
#'               'Wstrzymal sie','Nieobecny',
#'     6) club - deputy's club,
#' 4. statements with columns:
#'     1) id_statement - statement's id, like: 
#'     (meeting's number).(voting's number).(statement's number),
#'     2) nr_term_of_office - Polish Diet's number of term of office,
#'     3) surname_name - author of statement,
#'     4) date_statement - statement's date,
#'     5) titles_order_points - title of order points,
#'     6) statement - content of statement.}
#' 
#' @usage create_database(dbname, user, password, host)
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
#' create_database(dbname, user, password, host)}
#' 
#' @note
#' All information is stored in PostgreSQL database.
#'
#' @author Piotr Smuda
#' 
#' @export
#'
#' @import RPostgreSQL rvest stringi
#' @importFrom DBI dbDriver
#' @importFrom XML readHTMLTable
#' 

create_database <- function(dbname, user, password, host) {
    stopifnot(is.character(dbname), is.character(user), is.character(password), is.character(host))
    
    # connecting to database
    drv <- dbDriver("PostgreSQL")
    database_diet <- dbConnect(drv, dbname = dbname, user = user, password = password, host = host)
    
    # creating table with deputies data
    dbSendQuery(database_diet, "CREATE TABLE deputies (id_deputy varchar(4) NOT NULL,
                nr_term_of_office int NOT NULL, surname_name varchar(50) NOT NULL,
                PRIMARY KEY (id_deputy, nr_term_of_office),
                CONSTRAINT uq_surname_name UNIQUE (nr_term_of_office, surname_name))")
    
    # creating table with voting data
    dbSendQuery(database_diet, "CREATE TABLE votings (id_voting int NOT NULL,
                nr_term_of_office int NOT NULL, nr_meeting int NOT NULL,
                date_meeting date NOT NULL, nr_voting int NOT NULL,
                topic_voting text NOT NULL, link_results varchar(200),
                PRIMARY KEY (id_voting, nr_term_of_office))")
    
    # creating table with votes data
    dbSendQuery(database_diet, "CREATE TABLE votes (id_vote int NOT NULL, 
                nr_term_of_office int NOT NULL, id_deputy varchar(4) NOT NULL,
                id_voting int NOT NULL, vote varchar(20) NOT NULL, club varchar(50),
                PRIMARY KEY (id_vote, nr_term_of_office),
                FOREIGN KEY (id_deputy, nr_term_of_office) REFERENCES deputies(id_deputy, nr_term_of_office),
                FOREIGN KEY (id_voting, nr_term_of_office) REFERENCES votings(id_voting, nr_term_of_office))")
    
    # creating table with statements data
    dbSendQuery(database_diet, "CREATE TABLE statements (id_statement varchar(11) NOT NULL,
                nr_term_of_office int NOT NULL, surname_name varchar(100) NOT NULL,
                date_statement date NOT NULL, titles_order_points text NOT NULL,
                statement text NOT NULL,
                PRIMARY KEY (id_statement, nr_term_of_office))")
    
    # creating table with counter data
    dbSendQuery(database_diet, "CREATE TABLE counter (id SERIAL PRIMARY KEY, what varchar(10) NOT NULL,
                date varchar(10) NOT NULL)")
    
    # disconnecting to database
    suppressWarnings(dbDisconnect(database_diet))
    return(invisible(NULL))
} 

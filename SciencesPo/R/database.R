#' Wrapper for dbConnect
#'
#' Connects to a SQLite database or creates one if it does not already exist.
#'
#' If the '.sqlite' file extension is ommited from the dbname argument it is automatically added.
#'
#' @export
#'
#' @param dbname A character name or path to database file.
#' @return SQLiteConnection object.
#' @importFrom  stringr str_detect
#' @importFrom RSQLite dbConnect
#' @examples \dontrun{
#' db <- database("mydb")
#' }
database <- function(dbname){
  if(!stringr::str_detect(dbname, "\\.sqlite$")) {
    dbname <- paste(dbname, "sqlite", sep = ".")
  }
  RSQLite::dbConnect(RSQLite::SQLite(), dbname)
}
NULL


#' Creates a temporary table in the database
#'
#' This function is useful if most of your work is on a subset of the database
#'
#' The table will exist for as long as the database connection is kept open
#' The Select_query argument will take the output from a select_events(sql_only = TRUE) based function
#'
#' @export
#'
#' @param db a database connection object
#' @param tab_name character name for the teporary table
#' @param query character the query that specifies the temporary table
#' @examples \dontrun{
#' db <- database("Eleicoes")
#' dbTempTable(db, tab_name = "gerais_2006",
#'            query = select_events(db, tab = "Candidatos",
#'            columns = c("nome", "resultado", "partido"),
#'            where = "data > '1997-01-01'",
#'            sql_only = TRUE))
#' }
#' @importFrom  RSQLite dbGetQuery
#'
dbTempTable <- function(db, tab_name, query){
  RSQLite::dbGetQuery(db, paste("CREATE TEMP TABLE", tab_name, "AS", query, ";"))
  message(sprintf("Temporary table '%s' created", tab_name))
}



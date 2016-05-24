#' Conveniently Create Connection Object to PostgreSQL based timeseriesdb
#' 
#' Create a conection object while getting user information from the R session. 
#' Also standard db parameters like port and driver are set. Yet flexible information like 
#' host or dbname should be added to Sys.setenv environments. 
#' 
#' @param dbuser character username. Defaults to reading username from Sys.info()
#' @param dbname character name of the database, assumes dbname is stored in TIMESERIESDB_NAME.
#' @param dbhost character host address, asssumes dbhost ist stored in TIMESERIESDB_HOST.
#' @param passwd character password is used. No defaults, best way to pass a password is to 
#' .rs.askForPassword to hide password entries when using R Studio.
#' @param dbport integer port number defaults to 5432 for postgres
#' @export
createConObj <- function(dbuser = Sys.info()["user"],
                         dbname = Sys.getenv("TIMESERIESDB_NAME"),
                         dbhost = Sys.getenv("TIMESERIESDB_HOST"),
                         passwd,
                         dbport = 5432){
  drv = DBI::dbDriver("PostgreSQL")
  
  con <- DBI::dbConnect(drv,
                        host = dbhost,
                        dbname = dbname,
                        port = dbport,
                        user = dbuser,
                        password = passwd)
  con
}



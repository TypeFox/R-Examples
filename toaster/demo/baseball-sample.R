#' @demoTitle baseball-sample
#' 
#' Demo sampling data from Aster tables
#'
#' To install and use baseball demo dataset in Aster:
#'
#' 1. download baseball.zip from
#'   https://bitbucket.org/grigory/toaster/downloads/baseball.zip
#' 2. run script to create data set in Aster
#'   sh load_baseball_data.sh -d mydbname -U username -w mypassword 
#' 3. create Aster ODBC DSN on your desktop
#'   see https://bitbucket.org/grigory/toaster/wiki/Home#markdown-header-odbc-driver-and-dns

library(toaster)

## utility input function
readlineDef <- function(prompt, default) {
  if (!is.null(prompt))
    prompt = paste0(prompt, "[", default, "]: ")
  else 
    prompt = paste0(prompt, ": ")
  
  result = readline(prompt)
  if (result == "") 
    return (default)
  else
    return (result)
}

## utility connection function
connectWithDSNToAster <- function(dsn=NULL) {
  dsn = readlineDef("Enter Aster ODBC DSN: ", dsn)
  
  tryCatch(close(conn), error=function(err) {NULL})
  
  conn = tryCatch({
    conn = odbcConnect(dsn)
    odbcGetInfo(conn)
    return (conn)
  }, error=function(err) {
    stop(paste("Can't connect to Aster - check DSN '", dsn, "'"))
  })
}

## utility pause function
pause <- function() {
  cat("Press ENTER/RETURN/NEWLINE to continue.")
  readLines(n=1)
  invisible()
}

## connect to Aster first
conn = connectWithDSNToAster()

## must be connected to baseball dataset
if(!all(isTable(conn, c('pitching', 'batting')))) {
  stop("Must connect to baseball dataset and tables must exist.")
}

# Sample 1% from batting table
batters = computeSample(conn, "batting", sampleFraction=0.01)
dim(batters)
head(batters)

pause()

# Sample 1K of pitching rows for AL only
pitchers = computeSample(conn, "pitching", sampleSize=1000, where="lgid='AL'")
dim(pitchers)
tail(pitchers)

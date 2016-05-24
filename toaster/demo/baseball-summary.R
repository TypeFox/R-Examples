#' @demoTitle baseball-summary
#' 
#' Demo table summary statistics
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
if(!all(isTable(conn, c('pitching_enh', 'batting_enh', 'teams_enh', 'master_enh')))) {
  stop("Must connect to baseball dataset and tables must exist.")
}


## compute pitching summary statistics 
pitchingInfoDemo = getTableSummary(channel=conn, 'pitching_enh', 
                                   except=c('lgid', 'teamid', 'playerid', 'yearid', 'decadeid'))
## Summary stats computed for the columns:
pitchingInfoDemo$COLUMN_NAME
## Column attributes and summary statistics collected:
names(pitchingInfoDemo)

pause()

## Compute statistics on subset of data:
## batting in 2000s
battingInfo2000s = getTableSummary(channel=conn, 'batting_enh', 
                                   where='yearid between 2000 and 2012')

pause()
 
## compute statistics on certain pitching metrics 
## including every 3d percentile from 1 to 99 and the modes
pitchingInfoWithMode = getTableSummary(channel=conn, 'pitching_enh',
                             include=c('hr', 'bb', 'so', 'lgid'),
                             percentiles=seq(1,99,3), mode=TRUE)
## all computed statistics
names(pitchingInfoWithMode)


pause()
                              

## compute summary statitics on all numeric columns
## except few
teamInfoNumeric = getTableSummary(channel=conn, 'teams_enh', 
                           include=getNumericColumns(sqlColumns(conn, 'teams_enh')),
                           except=c('lgid', 'teamid', 'playerid', 'yearid', 'decadeid'))

pause()


## compute statistics on temporal columns
masterTemporalInfo = getTableSummary(channel=conn, 'master_enh', 
                                     include=getTemporalColumns(sqlColumns(conn, 'master_enh')))
## all computed statistics:
names(masterTemporalInfo)
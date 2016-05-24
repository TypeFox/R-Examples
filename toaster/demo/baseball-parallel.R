#' @demoTitle baseball-parallel
#' 
#' Demo parallel table summary and percentiles computations in R
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
if(!all(isTable(conn, c('batting', 'batting_enh')))) {
  stop("Must connect to baseball dataset and tables must exist.")
}

## compute table summary
system.time(
  battingInfo <- getTableSummary(channel=conn, 'batting',
                                 include=c('r', 'h', 'g', 'hr', 'bb', 'so', 'x2b', 'x3b', 
                                           'lgid', 'playerid', 'yearid', 'teamid'),
                                 mode=TRUE, parallel=FALSE)
)

pause()

## the same parallel 
library(doParallel)

cl = makeCluster(6)
registerDoParallel(cl, cores=6)
system.time(
  battingInfoPar <- getTableSummary(channel=conn, 'batting',
                                    include=c('r', 'h', 'g', 'hr', 'bb', 'so', 'x2b', 'x3b',
                                              'lgid', 'playerid', 'yearid', 'teamid'),
                                    mode=TRUE, parallel=TRUE)
)
stopCluster(cl)

pause()

## compute percentiles
system.time(
  allBattingPerc <- computePercentiles(conn, "batting_enh", columns=c("ba","g","r","h","x2b","x3b","hr","bb","so"),
                                       parallel=FALSE)
)

pause()

## the same parallel
cl = makeCluster(6)
registerDoParallel(cl, cores=6)
system.time(
  allBattingPerc <- computePercentiles(conn, "batting_enh", columns=c("ba","g","r","h","x2b","x3b","hr","bb","so"),
                                       parallel=TRUE)
)
stopCluster(cl)

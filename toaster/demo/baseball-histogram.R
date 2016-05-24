#' @demoTitle baseball-histogram
#' 
#' Demo histograms
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

## connect to Aster first
conn = connectWithDSNToAster()

## must be connected to baseball dataset
if(!all(isTable(conn, c('pitching_enh', 'batting_enh')))) {
  stop("Must connect to baseball dataset and tables must exist.")
}

# Compare Yankees and Rangers ERA in 2000s
h2000s = computeHistogram(conn, 'pitching_enh', 'era',
                          binsize=0.2, startvalue=0, endvalue=10, 
                          where="yearID between 2000 and 2013 and teamid in ('NYA','TEX')", by='teamid')
createHistogram(h2000s, fill='teamid', facet='teamid', title='TEX vs. NYY (ERA, 2000-2013)', xlab='ERA', ylab='count',
                legendPosition='none')


# Compare League BA by decades
hBADecadeLg=computeHistogram(conn, 'batting_enh', 'ba', by=c('lgid','decadeid'),
                             startvalue=0.01, endvalue=0.5, binsize=0.01,
                             where="decadeid between 1960 and 2000")
createHistogram(hBADecadeLg, fill='lgid', facet=c('decadeid','lgid'), facetScales='fixed', 
                title='AL vs. NL BA 1960-2009', xlab='BA', ylab='count',
                breaks=seq(from=0.0, to=0.5, by=0.05), paletteValues = c("red","blue"),
                legendPosition='none')

# Compare League ERA by decades
hERADecadeLg=computeHistogram(conn, 'pitching_enh', 'era',
                              startvalue=0.01, endvalue=10., binsize=0.1,
                              where="decadeid between 1960 and 2000", by=c('lgid','decadeid'))
createHistogram(hERADecadeLg, fill='lgid', facet=c('lgid','decadeid'), facetScales='fixed', 
                title='AL vs. NL ERA 1960-2009', xlab='ERA', ylab='count',
                breaks=seq(from=0.0, to=10.0, by=0.5), paletteValues = c("red","blue"),
                legendPosition='none')

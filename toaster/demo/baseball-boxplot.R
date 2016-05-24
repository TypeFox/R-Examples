#' @demoTitle baseball-boxplots
#' 
#' Demo percentiles and box plot visuals
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
if(!all(isTable(conn, c('master_enh', 'batting_enh')))) {
  stop("Must connect to baseball dataset and tables must exist.")
}

## Single boxplot: BA of all players
allBAPerc = computePercentiles(conn, "batting_enh", columns="ba")
createBoxplot(allBAPerc, fill=NULL, title="BA Boxplot")

## Team BA in 2000s 
teamBA2000sPerc = computePercentiles(conn, "batting_enh", columns="ba",
                                     by="teamid", where="yearid >= 2000")
createBoxplot(teamBA2000sPerc, x="teamid", fill="teamid", useIQR=TRUE,
              title="Team BA in 2000s", xlab="Team", ylab=NULL,
              coordFlip=TRUE, legendPosition="none")
## and without using IQR
createBoxplot(teamBA2000sPerc, x="teamid", fill="teamid", 
              title="Team BA in 2000s (no IQR)", xlab="Team", ylab=NULL,
              coordFlip=TRUE, legendPosition="none")

## Team BA by League facets in 2000s
teamBA2000byLgPerc = computePercentiles(conn, "batting_enh", columns="ba", 
                                        by=c('teamid','lgid'), where="yearid >= 2000")
createBoxplot(teamBA2000byLgPerc, x="teamid", fill="lgid", useIQR=TRUE,
              facet=c("lgid"), paletteValues = c("red","blue"), 
              title="Team BA by League in 2000s",
              coordFlip=TRUE, legendPosition="none")

## League BA by decade and league facet grid
teamsBAbyDecadePerc = computePercentiles(conn, "batting_enh", columns="ba",
                                         by=c("lgid","decadeid"), 
                                         where="yearid>=1960 and lgid in ('AL','NL')")
createBoxplot(teamsBAbyDecadePerc, facet=c("lgid","decadeid"), fill="lgid", useIQR=TRUE,
              paletteValues = c("red","blue"), legendPosition="none",
              title="League BA by Decade")

## League BA, SLG, TA, RBI
teamsAllBatting = computePercentiles(conn, "batting_enh", columns = c("ba", "slg", "ta"),
                              by=c('lgid','decadeid'), where="yearid >= 1980")
createBoxplot(teamsAllBatting, x='column', facet=c('lgid', 'decadeid'), useIQR=TRUE, coordFlip=TRUE, 
              title="Batting by League and Decades",
              legendPosition="none")

# Temporal percentiles (on epoch value)
playerAllDates = computePercentiles(conn, "master_enh", columns=c('debut','finalgame','birthdate','deathdate'),
                                    temporal=TRUE, percentiles=c(0))
createBoxplot(playerAllDates, x='column', value='epoch', useIQR=TRUE, 
              title="Boxplots for Date columns (epoch values)", legendPosition="none")

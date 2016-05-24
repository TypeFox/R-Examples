#' @demoTitle baseball-heatmap
#' 
#' Demo heatmaps
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
library(scales)

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
if(!all(isTable(conn, c('pitching_enh', 'batting_enh', 'teams_enh')))) {
  stop("Must connect to baseball dataset and tables must exist.")
}

# Heatmap of Average Wins by Franchise and Decade
hmFranchWins = computeHeatmap(conn, "teams_enh", 'franchid', 'decadeid', 'avg(w) w', 
                    where="decadeid >= 1950")
createHeatmap(hmFranchWins, 'decadeid', 'franchid', 'w',
              title = "Average Wins by Franchise and Decade")
                     
# Heatmap of Average Diff W-L with diverging color gradient
hmFranchWinLoss = computeHeatmap(conn, "teams_enh", 'franchid', 'decadeid', 'avg(w-l) wl', 
                    where="decadeid >= 1950")
createHeatmap(hmFranchWinLoss, 'decadeid', 'franchid', 'wl', divergingColourGradient = TRUE,
              title = "Average Win-Loss by Franchise and Decade")

# BA Heatmaps by Leagues
heatBA = computeHeatmap(conn, 'batting_enh', 'teamid', 'yearid',
                        aggregates="SUM(BA*AB)/SUM(AB) ba",
                        where="yearid between 1990 and 2012", by="lgid",
                        withMelt=FALSE)
createHeatmap(heatBA, 'yearid', 'teamid', 'ba',
              title='BA Heatmap AL vs. NL 1990-2012', xlab='Year', ylab='Team',
              lowGradient="lightgrey", highGradient="red",
              text=TRUE, percent=FALSE, digits=3,
              facet='lgid', legendPosition="none")

# ERA Heatmaps by Leagues
pitchERA = computeHeatmap(conn, 'pitching_enh', 'teamid', 'yearid',
                         aggregates=c("SUM(IPOuts*ERA)/SUM(IPOuts) era"),
                         where="yearid between 1990 and 2012", by="lgid")
createHeatmap(pitchERA, 'yearid', 'teamid', 'era',
              title='ERA Heatmap AL vs. NL 1990-2012', xlab='Year', ylab='Team',
              lowGradient="red", highGradient="lightgrey",
              text=TRUE, percent=FALSE, digits=3,
              facet='lgid', legendPosition="none")

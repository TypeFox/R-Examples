#' @demoTitle baseball-barchart
#' 
#' Demo bar charts
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
library(ggplot2)

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
if(!all(isTable(conn, c('pitching_enh', 'teams_enh')))) {
  stop("Must connect to baseball dataset and tables must exist.")
}

# simple bar chart
bc = computeBarchart(channel=conn, tableName="pitching_enh", category="teamid",
                     aggregates="AVG(era) era", 
                     where="yearid >= 2000 and lgid='AL'")
bc = transform(bc, teamid2=reorder(factor(teamid), rank(era))) # reorder by era
createHistogram(bc, "teamid2", "era", fill="teamid2", 
                title = "AL Teams Average ERA in 2000s", legendPosition="none")

# multipe aggregates in the same bar chart (with melt)
bc = computeBarchart(channel=conn, tableName="pitching_enh", category="teamid",
                     aggregates=c("AVG(era) era", "AVG(whip) whip"), withMelt=TRUE,
                     where="yearid >= 2000 and lgid='AL'")
bc = transform(bc, teamid2=reorder(factor(teamid), rank(value))) # reorder by era
createHistogram(bc, "teamid2", "value", fill="teamid2", facet="variable",
                title = "AL Teams Average ERA and WHIP in 2000s", legendPosition="none")

# Average pitcher stats by team and decade
bc = computeBarchart(conn, "pitching_enh", "teamid", 
                     aggregates=c("AVG(era) era", "AVG(whip) whip", "AVG(ktobb) ktobb"),
                     where="yearid >= 1990 and lgid='AL'", by="decadeid", withMelt=TRUE)
bc = transform(bc, teamid2=reorder(factor(teamid), rank(value))) # reorder by value
createHistogram(bc, "teamid2", "value", fill="teamid2", facet=c("variable", "decadeid"), 
                legendPosition="bottom",
                title = "AL Teams Pitching Stats by decades (1990-2012)",
                themeExtra = guides(fill=guide_legend(nrow=2)))

# Franchise wins-loss Historical Trend by decades
franchwl = computeBarchart(conn, "teams_enh", "franchid",
                           aggregates=c("AVG(w-l) wl"),
                           by="decadeid",
                           where="yearid >=1960 and lgid = 'AL'")
createHistogram(franchwl, "decadeid", "wl", fillColour="white",
                facet="franchid", ncol=5, facetScales="free_y",
                legendPosition="none",
                trend=TRUE, trendLinesize=4,
                title="Franchise W-L Historical Trend by Decade (AL)",
                ylab="Average W-L")

# Players with top KtoBB average and average number of ipouts
playerTopKtoBB = computeBarchart(conn, "pitching_enh", "playerid",
               aggregates=c("AVG(ktobb) ktobb", "AVG(ipouts) ipouts"),
               where="lgid in ('AL','NL') and yearid > 1960",
               top=30, orderBy="ktobb desc")
createHistogram(playerTopKtoBB, "playerid", "ktobb", fill="ipouts",
                scaleGradient=scale_fill_gradient("IP Outs", high="red4", low="tan"),
                title="KtoBB by Team with IPOuts as Gradient Fill")

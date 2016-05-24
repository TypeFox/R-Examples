#' @demoTitle baseball-bubble
#' 
#' Demo bubble charts
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
if(!all(isTable(conn, c('teams_enh')))) {
  stop("Must connect to baseball dataset and table 'teams_enh' must exist.")
}

# Bubble chart example
bubble = computeHeatmap(conn, 'teams_enh', 'franchid', 'decadeid',
                        aggregates=c("SUM(BA*AB)/SUM(AB)ba", 
                                       "SUM(IPOuts*ERA)/SUM(IPOuts) era",
                                       "ROUND(AVG(8-rank)) rank",
                                       "MIN(lgid) lgid"),
                        where="yearid between 1970 and 2009")


createBubblechart(bubble, "ba", "era", "rank", label="franchid", fill="franchid",
                  facet=c("decadeid","lgid"), #ncol=1, 
                  scaleSize = FALSE, shapeSizeRange=c(1,15), shapeMaxSize = 15,
                  title="Team Ranks by BA and ERA", 
                  labelSize = 5, labelColour = "black", labelVJust = 1,
                  legendPosition="none", themeExtra = guides(fill = "legend", size = "legend"))

createBubblechart(bubble, "ba", "era", "rank", label="franchid", fill="franchid",
                  facet=c("lgid", "decadeid"), ncol=1, shapeSizeRange=c(1,20),
                  title="Team Ranks by BA and ERA", 
                  labelSize = 5, labelColour = "black", labelVJust = 1,
                  legendPosition="none", themeExtra = guides(fill = "legend", size = "legend"))
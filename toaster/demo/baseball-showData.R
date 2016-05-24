#' @demoTitle baseball-showData
#' 
#' Demo showData: 'overview', 'boxplot', 'corr' (correlation table), and 'scatterplot'
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

# utility test function
isTableInfo <- function(name) {
  return (exists(name) && class(get(name)) == "data.frame" && 
            all(c("TABLE_CAT","TABLE_SCHEM","TABLE_NAME","COLUMN_NAME","COLUMN_NAME", 
                  "DATA_TYPE","TYPE_NAME","COLUMN_SIZE","BUFFER_LENGTH","DECIMAL_DIGITS",
                  "NUM_PREC_RADIX","NULLABLE","REMARKS","COLUMN_DEF","SQL_DATA_TYPE",
                  "SQL_DATETIME_SUB","CHAR_OCTET_LENGTH","ORDINAL_POSITION","IS_NULLABLE",
                  "total_count","distinct_count","not_null_count","null_count","minimum",
                  "maximum","average","deviation","25%","50%","75%","IQR") %in% names(get(name))))
}

## connect to Aster first
conn = connectWithDSNToAster()

## must be connected to baseball dataset
if(!all(isTable(conn, c('pitching_enh', 'batting_enh', 'teams_enh')))) {
  stop("Must connect to baseball dataset and tables must exist.")
}

## compute table statistics as necessary
if (!isTableInfo("battingInfo"))
  pitchingInfo = getTableSummary(conn, 'pitching_enh', percentiles=c(10,25,50,75,90),
                                 where='yearid between 2000 and 2013')

if (!isTableInfo("pitchingInfo"))
  battingInfo = getTableSummary(conn, 'batting_enh', where='yearid between 2000 and 2013')

if (!isTableInfo("teamsInfo"))
  teamsInfo = getTableSummary(conn, tableName='teams_enh', where='yearid between 1960 and 2013')

## overview of pitching metrics,
## with free y-axis scales
showData(channel=NULL, 'pitching_enh', pitchingInfo, type='numeric', format='overview',
         include=c('h','so','r','g','bb','er','hr'), 
         measures = c('average','deviation','IQR','minimum','10%','25%','50%','75%','90%','maximum'),
         scales="free_y", where='yearid between 2000 and 2013',
         title="Overview format: Data Statistics by rows, Data Columns (Pitching Stats) by columns")


## box plots
showData(channel=NULL, tableName='pitching_enh', tableInfo=pitchingInfo, format='boxplot',
         except = c('yearid', 'decadeid'), coordFlip=TRUE,
         title="Boxplots of Pitching Stats")

showData(channel=NULL, tableName='pitching_enh', tableInfo=pitchingInfo, format='boxplot',
         include=c('bfp','er','h','ipouts','r','so'), ncol=3,
         facet=TRUE, scale="free_x", defaultTheme=theme_grey(),
         title="Pitching Stats Boxplots in Facets with Grey Theme")

## histograms
showData(conn, 'pitching_enh', tableInfo=pitchingInfo, format='histogram', 
         include=c('baopp', 'era', 'whip', 'ktobb', 'fip'), 
         numBins=50, facet=TRUE, where='yearid between 2000 and 2013',
         title="Histograms of Pitching Stats in Facets")

showData(conn, 'batting_enh', battingInfo,
         include=c('ba','ta','slg'),
         format='histogram', numBins=100, facet=TRUE,
         where='yearid between 2000 and 2013 
                and ab >= 30',
         title="Histograms of Batting Stats in 2000s with At Least 30 AB")

## correlation
showData(conn, tableName='pitching_enh', tableInfo=pitchingInfo, format='corr', 
         include=c('w', 'l', 'g', 'ipouts', 'h', 'er', 'hr', 'bb', 'so', 'bfp', 'r'),
         corrLabel='value', digits=2, shapeSizeRange=c(5,25), 
         defaultTheme=theme_classic(base_size=12), legendPosition="none",
         title="Correlation Matrix of Pitching Stats")

showData(conn, tableName='teams_enh', tableInfo=teamsInfo, format='corr', 
         except=c('yearid','decadeid','ghome','g'),
         corrLabel='value', digits=2, shapeSizeRange=c(5,25),
         defaultTheme=theme_classic(base_size=12), legendPosition="none",
         title="correlation Matrix of Team Stats")

## scatterplots
showData(conn, 'pitching_enh', tableInfo=pitchingInfo, format='scatterplot', 
         include=c('so', 'er'), pointColour="lgid", 
         sampleSize=10000, regressionLine=TRUE,
         where='yearid between 1980 and 2000', legendPosition="none",
         title="SO vs ER by League (1980-2000, AL: red, NL: grey)")

## the same with facets
showData(conn, 'pitching_enh', tableInfo=pitchingInfo, format='scatterplot', 
         include=c('so', 'er'), facetName="lgid", pointColour="lgid", 
         sampleSize=10000, regressionLine=TRUE,
         where='yearid between 1980 and 2000', legendPosition="none",
         title="SO vs ER by League (1980-2000)")

showData(conn, 'pitching_enh', format='scatterplot', tableInfo=pitchingInfo, 
         include=c('so','er'), facetName=c('lgid','decadeid'), pointColour="lgid",
         sampleFraction=0.25, regressionLine=TRUE,
         where='yearid between 1980 and 2013', legendPosition="none",
         title="SO vs ER by League and Decade (1980 - 2013)")

showData(conn, 'pitching_enh', format='scatterplot', tableInfo=pitchingInfo, 
         include=c('ktobb', 'fip', 'teamid','yearid'),
         sampleFraction=0.25, facetName=c('teamid','decadeid'), regressionLine=TRUE,
         where="yearid between 1970 and 2013 and teamid in ('TEX','NYA')",
         title="Pitching KTOBB vs FIP by Decade and Team (Yankees and Rangers)")


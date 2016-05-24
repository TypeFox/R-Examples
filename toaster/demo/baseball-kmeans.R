#' @demoTitle baseball-kmeans
#' 
#' Demo kmeans clustering of baseball batting Aster table
#'
#' To install and use baseball demo dataset in Aster:
#'
#' 1. download baseball.zip from
#'   https://raw.githubusercontent.com/wiki/teradata-aster-field/toaster/downloads/baseball.zip
#'   and extract files into directory of your choice
#' 2. in the directory run the script to create data set in Aster
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
if(!all(isTable(conn, c('batting')))) {
  stop("Must connect to baseball dataset and tables must exist.")
}

# run kmeans in Aster
km.demo = computeKmeans(conn, "batting_enh", centers=3, iterMax = 25,
                   aggregates = c("COUNT(*) cnt", "AVG(g) avg_g", "AVG(r) avg_r", "AVG(h) avg_h","AVG(ab) avg_ab", 
                                  "AVG(ba) ba", "AVG(slg) slg", "AVG(ta) ta"),
                   id="playerid || '-' || teamid || '-' || yearid", include=c('g','r','h','ab'),
                   # scaledTableName='kmeans_demo_scaled', centroidTableName='kmeans_demo_centroids', schema='public',
                   where="yearid > 2000", test=FALSE)

createCentroidPlot(km.demo, format="line")
createCentroidPlot(km.demo, format="line", groupByCluster=FALSE)

createCentroidPlot(km.demo, format="bar")
createCentroidPlot(km.demo, format="bar", groupByCluster=FALSE)

createCentroidPlot(km.demo, format="bar_dodge")
createCentroidPlot(km.demo, format="bar_dodge", groupByCluster=FALSE)

createCentroidPlot(km.demo, format="heatmap")
createCentroidPlot(km.demo, format="heatmap", coordFlip = TRUE)

createClusterPlot(km.demo)
createClusterPlot(km.demo, colorByCluster = FALSE)

pause()

# sample clustered data
km.demo = computeClusterSample(conn, km.demo, '0.5', test=FALSE)

createClusterPairsPlot(km.demo)

pause()


# silhouette analysis
km.demo = computeSilhouette(conn, km.demo)

createSilhouetteProfile(km.demo)

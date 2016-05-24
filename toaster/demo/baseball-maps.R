#' @demoTitle baseball-maps
#' 
#' Demo geo-spatial data and maps
#'
#' To install and use baseball demo dataset in Aster:
#'
#' 1. download baseball.zip from
#'   https://bitbucket.org/grigory/toaster/downloads/baseball.zip
#' 2. run script to create data set in Aster
#'   sh load_baseball_data.sh -d mydbname -U username -w mypassword 
#' 3. create Aster ODBC DSN on your desktop
#'   see https://bitbucket.org/grigory/toaster/wiki/Home#markdown-header-odbc-driver-and-dns

require(toaster)
require(memoise)
require(ggmap)

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

## utility function to test Internet connectivity
havingIP <- function() {
  if (.Platform$OS.type == "windows") {
    ipmessage <- system("ipconfig", intern = TRUE)
  } else {
    ipmessage <- system("ifconfig", intern = TRUE)
  }
  validIP <- "((25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)[.]){3}(25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)"
  any(grep(validIP, ipmessage))
}

## must have Internet
if(!havingIP()) {
  stop("Must have Internet connection to demo map functionality!")
}

## connect to Aster first
conn = connectWithDSNToAster()

## must be connected to baseball dataset
if(!all(isTable(conn, "teams_enh"))) {
  stop("Must connect to baseball dataset and table 'teams_enh' must exist.")
}

data = computeAggregates(conn, "teams_enh",
               by = c("name || ', ' || park teamname", "lgid", "teamid", "decadeid"),
               aggregates = c("min(name) name", "min(park) park", "10-avg(rank) rank", 
                              "avg(attendance) attendance"))

if (!exists("geocodeMem")) {
  geocodeMem = memoise(geocode) 
}

createMap(data=data[data$decadeid>=2000,], source = "google", maptype = "roadmap", zoom=4, 
              facet=c("lgid", "decadeid"),
              locationName='teamname', metricName='attendance', labelName='name',
              shapeColour="blue", scaleRange = c(2,12), textColour="black",
              title='Team Attendance by Decade and Leage',
              geocodeFun=geocodeMem)

createMap(data=data[data$decadeid==2000 & data$lgid == 'AL',], source = "stamen", maptype = "toner", zoom = 4, 
          locationName='teamname', metricName='rank', labelName='name',
          shapeColour="blue", scaleRange=c(4,20), textColour="black",
          title='Average AL Team Rank in 1990s (higher is better)',
          geocodeFun=geocodeMem)
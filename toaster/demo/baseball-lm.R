#' @demoTitle baseball-lm
#' 
#' Demo linear regression model (lm)
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

## utility connection functions
connectAdHocToAster <- function(uid="beehive", pwd=NULL, server="localhost", port="2406", database="beehive") {
  uid = readlineDef("Enter user id: ", uid)
  pwd = readlineDef("Enter password: ", pwd)
  server = readlineDef("Enter server: ", server)
  port = readlineDef("Enter port: ", port)
  database = readlineDef("Enter database: ", database)
  
  tryCatch(close(conn), error=function(err) {NULL})
  
  connectStr = paste0("driver={Aster ODBC Driver};server=", server, ";port=", port, 
                      ";database=",database,";uid=",uid,";pwd=", pwd)
  conn = tryCatch({
    
    conn = odbcDriverConnect(connection=connectStr)
    odbcGetInfo(conn)
    return (conn)
  }, error=function(err) {
    stop("Can't connect to Aster - check ip/port/database/uid/pwd", err)
  })
}

connectWithDSNToAster <- function(dsn=NULL) {
  dsn = readlineDef("Enter Aster ODBC DSN: ", dsn)
  
  tryCatch(close(conn), error=function(err) {NULL})
  
  conn = tryCatch({
    conn = odbcConnect(dsn)
    odbcGetInfo(conn)
    return (conn)
  }, error=function(err) {
    stop(paste("Can't connect to Aster - check DSN:", dsn))
  })
}

pause <- function() {
  cat("Press ENTER/RETURN/NEWLINE to continue.")
  readLines(n=1)
  invisible()
}

### connect first
conn = connectWithDSNToAster()

## must be connected to baseball dataset
if(!all(isTable(conn, c('pitching_enh', 'batting_enh')))) {
  stop("Must connect to baseball dataset and tables must exist.")
}

### simple model with 3 numerical predictors
model1 = computeLm(channel=conn, tableName="batting_enh",
                   formula= ba ~ rbi + bb + so, sampleSize=10000)
summary(model1)
plot(model1)

### compare with lm()
data = computeSample(channel=conn, tableName="batting_enh",
                     include=c("ba","rbi","bb","so"), sampleSize=10000)
fit1 = lm(formula= ba ~ rbi + bb + so, data=data)
summary(fit1)
plot(fit1)

pause()

### added ER predictor to the model and defined subset with where
modelNL = computeLm(channel=conn, tableName="pitching_enh", formula= era ~ er + hr + bb + so, 
                    where = "yearid >= 2000 and lgid = 'NL'")
summary(modelNL)
plot(modelNL)

pause()

### added categorical predictor
modelLg = computeLm(channel=conn, tableName="batting_enh", formula=ba ~ rbi + bb + so + lgid, 
                    where="lgid in ('AL','NL')")
summary(modelLg)

pause()

### category with more than 2 values, also lowered sample size to 500 rows 
### (coefficients are always computed on all data)
modelTeam10K = computeLm(channel=conn, tableName="batting_enh", formula=ba ~ rbi + bb + so + teamid,
                    sampleSize = 50, where="teamid in ('TEX','NYY','OAK','PIT','DET') and yearid >= 1990")
summary(modelTeam10K)

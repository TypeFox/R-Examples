service <- Sys.getenv("_R_CHECK_HAVE_MYSQL_")
require("tfplot")

if(identical(as.logical(service), TRUE)) {

require("TSMySQL")
require("timeSeries")

cat("***** RMySQL  with timeSeries representation *********\n")
dbname   <- Sys.getenv("MYSQL_DATABASE")
if ("" == dbname)   dbname <- "test"

user	<- Sys.getenv("MYSQL_USER")

if ("" != user) {
    # specifying host as NULL or "localhost" results in a socket connection
    host    <- Sys.getenv("MYSQL_HOST")
    if ("" == host)	host <- Sys.info()["nodename"] 
    passwd  <- Sys.getenv("MYSQL_PASSWD")
    if ("" == passwd)	passwd <- NULL
    #  See  ?"dbConnect-methods"
    setup <- RMySQL::dbConnect(RMySQL::MySQL(),
       username=user, password=passwd, host=host, dbname=dbname)  
  }else setup <- RMySQL::dbConnect(RMySQL::MySQL(), dbname=dbname) #user/passwd/host in ~/.my.cnf


TSsql::removeTSdbTables(setup, yesIknowWhatIamDoing=TRUE)
TSsql::createTSdbTables(setup, index=FALSE)

  
if ("" != user)  
  con <- tryCatch(TSconnect("MySQL", dbname=dbname, username=user, password=passwd, host=host))
else
  con <- tryCatch(TSconnect("MySQL", dbname=dbname)) # pass user/passwd/host in ~/.my.cnf

if(inherits(con, "try-error")) stop("Cannot connect to TS MySQL database.")

z <- ts(matrix(rnorm(10),10,1), start=c(1990,1), frequency=1)
TSput(z, serIDs="Series 1", con) 

# timeSeries seems to have a bug
#z2 <- timeSeries:::as.timeSeries.ts(z)
#z2 <- timeSeries:::as.timeSeries.zoo(zoo(z))
#z2 <- timeSeries:::as.timeSeries.zoo(zooreg(z))
#start(z2)
#time(z2)

# these retrieve but lose date info
z <- TSget("Series 1", con, TSrepresentation="timeSeries")
if("timeSeries" != class(z)) stop("timeSeries class object not returned.")

z <- TSget("Series 1", con, TSrepresentation=timeSeries::timeSeries)
if("timeSeries" != class(z)) stop("timeSeries class object not returned.")

#unresolved problem with timeSeries  (see above)
#tfplot(z)
#gives  'origin' must be supplied

TSrefperiod(z) 
TSdescription(z) 

#unresolved problem 
#tfplot(z, start="1991-01-01", Title="Test")
#gives  'origin' must be supplied

cat("**************        remove test tables\n")
TSsql::removeTSdbTables(con, yesIknowWhatIamDoing=TRUE)

cat("**************        disconnecting test\n")
dbDisconnect(con)

} else  {
   cat("MYSQL not available. Skipping tests.\n")
   cat("_R_CHECK_HAVE_MYSQL_ setting ", service, "\n")
   }

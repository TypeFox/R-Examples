# sourcing HistQuote.TSsql requires the Internet

service <- Sys.getenv("_R_CHECK_HAVE_ODBC_")

if(identical(as.logical(service), TRUE)) {

cat("************** RODBC  Examples ******************************\n")
cat("**************************************************************\n")
cat("* WARNING: THIS OVERWRITES TABLES IN TEST DATABASE ON SERVER**\n")
cat("**************************************************************\n")

###### This is to set up tables. Otherwise use TSconnect#########
  # This will fail if ODBC support is not installed or set up on the system 
  # with a message like:  [RODBC] ERROR: state IM002, code 0, 
  # message [unixODBC][Driver Manager]Data source name not found, 
  # and no default driver specified

   dbname   <- Sys.getenv("ODBC_DATABASE")
   if ("" == dbname)   dbname <- "test"

   require("RODBC")
   cat("**********setup 0\n")
   user    <- Sys.getenv("ODBC_USER")
   if ("" != user) {
       host    <- Sys.getenv("ODBC_HOST")
       if ("" == host)     host <- Sys.info()["nodename"] 
       passwd  <- Sys.getenv("ODBC_PASSWD")
       if ("" == passwd)   passwd <- NULL
       #  See  ?odbcConnect   ?odbcDriverConnect
       setup <- RODBC::odbcConnect(dsn=dbname, uid=user, pwd=passwd) #, connection=host) 
     }else  
       setup <- RODBC::odbcConnect(dsn=dbname) # pass user/passwd/host in ~/.odbc.ini

require("TSodbc") # because next needs DBI faking method dbExistsTable, etc

cat("**********setup 1\n")
# ToLower=TRUE here because database is PostgreSQL
TSsql::removeTSdbTables(setup, yesIknowWhatIamDoing=TRUE, ToLower=TRUE)

cat("**********setup 2\n")
TSsql::createTSdbTables(setup, index=FALSE)

cat("**********setup 3\n")
RODBC::odbcClose(setup)

detach(package:RODBC)

##################################################################
require("TSodbc")

con <- if ("" != user)  
          tryCatch(TSconnect("odbc", dbname=dbname, username=user, password=passwd, host=host)) 
    else  tryCatch(TSconnect("odbc", dbname=dbname)) # pass user/passwd/host in config file

if(inherits(con, "try-error")) cat("CreateTables did not work.\n")
 else {
   source(system.file("TSsql/Populate.TSsql", package = "TSsql"))
   source(system.file("TSsql/TSdbi.TSsql", package = "TSsql"))
   source(system.file("TSsql/dbGetQuery.TSsql", package = "TSsql"))
   source(system.file("TSsql/HistQuote.TSsql", package = "TSsql"))

   cat("**************        removing test database tables\n")
   TSsql::removeTSdbTables(con, yesIknowWhatIamDoing=TRUE, ToLower=TRUE)

   cat("**************        disconnecting test\n")
   dbDisconnect(con)
   }
} else {
   cat("ODBC not available. Skipping tests.\n")
   cat("_R_CHECK_HAVE_ODBC_ setting ", service, "\n")
   }

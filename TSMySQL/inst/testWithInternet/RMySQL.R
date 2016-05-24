# sourcing HistQuote.TSsql requires the Internet

service <- Sys.getenv("_R_CHECK_HAVE_MYSQL_")

if(identical(as.logical(service), TRUE)) {

cat("************** RMySQL  Examples ******************************\n")
cat("**************************************************************\n")
cat("* WARNING: THIS OVERWRITES TABLES IN TEST DATABASE ON SERVER**\n")
cat("**************************************************************\n")

###### This is to set up tables. Otherwise use TSconnect#########
   dbname   <- Sys.getenv("MYSQL_DATABASE")
   if ("" == dbname)   dbname <- "test"

   user    <- Sys.getenv("MYSQL_USER")
   if ("" != user) {
       # specifying host as NULL or "localhost" results in a socket connection
       host    <- Sys.getenv("MYSQL_HOST")
       if ("" == host)     host <- Sys.info()["nodename"] 
       passwd  <- Sys.getenv("MYSQL_PASSWD")
       if ("" == passwd)   passwd <- NULL
       #  See  ?"dbConnect-methods"
       setup <- RMySQL::dbConnect(RMySQL::MySQL(),
           username=user, password=passwd, host=host, dbname=dbname)  
     }else setup <- RMySQL::dbConnect(RMySQL::MySQL(), dbname=dbname) #user/passwd/host in ~/.my.cnf

DBI::dbListTables(setup) 

TSsql::removeTSdbTables(setup, yesIknowWhatIamDoing=TRUE)
TSsql::createTSdbTables(setup, index=FALSE)

DBI::dbListTables(setup) 
DBI::dbDisconnect(setup)
##################################################################
require("TSMySQL")

con <- if ("" != user)  
          tryCatch(TSconnect("MySQL", dbname=dbname, username=user, password=passwd, host=host)) 
    else  tryCatch(TSconnect("MySQL", dbname=dbname)) # pass user/passwd/host in ~/.my.cnf

if(inherits(con, "try-error")) stop("CreateTables did not work.")

source(system.file("TSsql/Populate.TSsql", package = "TSsql"))
require("zoo")
require("tframePlus")
source(system.file("TSsql/TSdbi.TSsql", package = "TSsql"))

source(system.file("TSsql/dbGetQuery.TSsql", package = "TSsql"))
source(system.file("TSsql/HistQuote.TSsql", package = "TSsql"))

cat("**************        remove test tables\n")
TSsql::removeTSdbTables(con, yesIknowWhatIamDoing=TRUE)

cat("**************        disconnecting test\n")
dbDisconnect(con)

} else  {
   cat("MYSQL not available. Skipping tests.\n")
   cat("_R_CHECK_HAVE_MYSQL_ setting ", service, "\n")
   }

# sourcing HistQuote.TSsql requires the Internet

service <- Sys.getenv("_R_CHECK_HAVE_POSTGRES_")

if(identical(as.logical(service), TRUE)) {


cat("************** RPostgreSQL  Examples ******************************\n")
cat("**************************************************************\n")
cat("* WARNING: THIS OVERWRITES TABLES IN TEST DATABASE ON SERVER**\n")
cat("**************************************************************\n")

###### This is to set up tables. Otherwise use TSconnect#########

   dbname   <- Sys.getenv("POSTGRES_DATABASE")
   if ("" == dbname)   dbname <- "test"

   user    <- Sys.getenv("POSTGRES_USER")
   host <- Sys.getenv("POSTGRES_HOST")
   if ("" == host) host  <- Sys.getenv("PGHOST")
   if ("" == host) host  <- "localhost"  #Sys.info()["nodename"] 
   if ("" != user) {
       passwd  <- Sys.getenv("POSTGRES_PASSWD")
       #  See  ?"dbConnect-methods"
       setup <- RPostgreSQL::dbConnect(RPostgreSQL::PostgreSQL(), dbname=dbname,
          user=user, password=passwd, host=host)  
     }else  {
	#( the postgres driver may also use PGDATABASE, PGHOST, PGPORT, PGUSER )
       # The Postgress documentation seems to suggest that it should be
       #   possible to get the host from the .pgpass file too, but I cannot.
       #get user/passwd in ~/.pgpass
       setup <- RPostgreSQL::dbConnect(RPostgreSQL::PostgreSQL(), dbname=dbname, host=host) 
       }

DBI::dbListTables(setup) 

TSsql::removeTSdbTables(setup, yesIknowWhatIamDoing=TRUE, ToLower=TRUE)
TSsql::createTSdbTables(setup, index=FALSE)

DBI::dbListTables(setup) 
DBI::dbDisconnect(setup)
##################################################################

require("TSPostgreSQL")

# pass user/passwd in ~/.pgpass (but host defaults to PGHOST or localhost).

if ("" != user)  
      con <- tryCatch(TSconnect("PostgreSQL", dbname=dbname, 
                              user=user, password=passwd, host=host)) 
else  con <- tryCatch(TSconnect("PostgreSQL", dbname=dbname)) 
    
if(inherits(con, "try-error")) stop("CreateTables did not work.")

source(system.file("TSsql/Populate.TSsql", package = "TSsql"))
source(system.file("TSsql/TSdbi.TSsql", package = "TSsql"))

source(system.file("TSsql/dbGetQuery.TSsql", package = "TSsql"))
source(system.file("TSsql/HistQuote.TSsql", package = "TSsql"))

cat("**************        disconnecting test\n")
TSsql::removeTSdbTables(con, yesIknowWhatIamDoing=TRUE, ToLower=TRUE)
dbDisconnect(con)

} else  {
   cat("POSTGRES not available. Skipping tests.\n")
   cat("_R_CHECK_HAVE_POSTGRES_ setting ", service, "\n")
   }

# This just uses dbConnect (not TSconnect) to see if things work

# Note that host will default to localhost if neither environment variables
#   PGHOST or POSTGRES_HOST are set. (Only one of them is needed.)

service <- Sys.getenv("_R_CHECK_HAVE_POSTGRES_")

Sys.info()

if(!identical(as.logical(service), TRUE)) {
   cat("POSTGRES not available. Skipping tests.\n")
   cat("_R_CHECK_HAVE_POSTGRES_ setting ", service, "\n")
 }else {
   require("TSPostgreSQL") 
   m <- "PostgreSQL"

   dbname   <- Sys.getenv("POSTGRES_DATABASE")
   if ("" == dbname)  { 
      dbname <- "test" #PostgreSQL default is template1
      cat("dbname set to test.\n")
      }
   else cat("dbname set to:", dbname, " by env variable POSTGRES_DATABASE\n") 

   host    <- Sys.getenv("POSTGRES_HOST")
   cat("host set to:", host, " by env variable POSTGRES_HOST\n") 
   
   user    <- Sys.getenv("POSTGRES_USER")
   if ("" != user) {
       cat("user set to:", user, " by env variable POSTGRES_USER\n") 
       # specifying host as NULL or "localhost" results in a socket connection
       if ("" == host) {
          host <- Sys.info()["nodename"]
          cat("host reset to nodename:", host, "\n") 
	  } 
       passwd  <- Sys.getenv("POSTGRES_PASSWD")
       cat("passwd set by env variable POSTGRES_PASSWD\n") 

       #  See  ?"dbConnect-methods"
       con <- RPostgreSQL::dbConnect(m, dbname=dbname,
          user=user, password=passwd, host=host)  
     }else  {
       cat("user set to empty string.\n") 
       #(postgres driver may also use PGDATABASE, PGHOST, PGPORT, PGUSER )
       # The Postgress documentation seems to suggest that it should be
       #   possible to get the host from the .pgpass file too, but I cannot.
       
       if ("" == host) {
          host <- Sys.getenv("PGHOST")
          cat("host reset to:", host, "by env variable PGHOST\n") 
	  } 
       if ("" == host) {
          host <- "localhost"  #Sys.info()["nodename"]
          cat("host reset to:", host, "\n") 
	  } 
       #get user/passwd in ~/.pgpass
       con <- RPostgreSQL::dbConnect(m, dbname=dbname, host=host) 
       }
   # dbListTables(con) needs a non-null dbname
   RPostgreSQL::dbDisconnect(con)
 }

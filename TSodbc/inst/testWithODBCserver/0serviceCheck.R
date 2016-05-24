# This just does a dbConnect (not a TSconnect) to see if things work

service <- Sys.getenv("_R_CHECK_HAVE_ODBC_")

Sys.info()

if(identical(as.logical(service), TRUE)) {
   #require("TSodbc") 
   require("RODBC") 

   dbname   <- Sys.getenv("ODBC_DATABASE")
   if ("" == dbname)   dbname <- NULL

   user    <- Sys.getenv("ODBC_USER")
   if ("" != user) {
       host    <- Sys.getenv("ODBC_HOST")
       if ("" != host) warning(
	 "host name is not used. The connection is taken from the .odbc.ini file")
       passwd  <- Sys.getenv("ODBC_PASSWD")
       if ("" == passwd)   passwd <- NULL
       #  See  ?odbcConnect   ?odbcDriverConnect
       con <- odbcConnect(dsn=dbname, uid=user, pwd=passwd) #, connection=host)  
    } else {
       if (is.null(dbname))   dbname <- "test"  # NULL dbname causes problem
       con <-  odbcConnect(dsn=dbname) # pass user/passwd/host in ~/.odbc.ini
       }
   odbcClose(con)
 } else {
   cat("ODBC not available. Skipping tests.\n")
   cat("_R_CHECK_HAVE_ODBC_ setting ", service, "\n")
   }

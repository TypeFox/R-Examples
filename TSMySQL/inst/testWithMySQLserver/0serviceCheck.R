# This just does a dbConnect (not a TSconnect) to see if things work

service <- Sys.getenv("_R_CHECK_HAVE_MYSQL_")

Sys.info()

if(identical(as.logical(service), TRUE)) {
   require("TSMySQL") 
   m <- RMySQL::MySQL()

   dbname   <- Sys.getenv("MYSQL_DATABASE")
   if ("" == dbname)   dbname <-  NULL 

   user    <- Sys.getenv("MYSQL_USER")
   if ("" != user) {
       # specifying host as NULL or "localhost" results in a socket connection
       host    <- Sys.getenv("MYSQL_HOST")
       if ("" == host)     host <- Sys.info()["nodename"] 
       passwd  <- Sys.getenv("MYSQL_PASSWD")
       if ("" == passwd)   passwd <- NULL
       #  See  ?"dbConnect-methods"
       con <- RMySQL::dbConnect(m,
          username=user, password=passwd, host=host, dbname=dbname)  
      }else {
        if (is.null(dbname))   dbname <- "test"  # NULL dbname causes segfault
	con <- RMySQL::dbConnect(m, dbname=dbname) # pass user/passwd/host in ~/.my.cnf
	}
   # dbListTables(con) needs a non-null dbname
   RMySQL::dbDisconnect(con)
 }else {
   cat("MYSQL not available. Skipping tests.\n")
   cat("_R_CHECK_HAVE_MYSQL_ setting ", service, "\n")
   }

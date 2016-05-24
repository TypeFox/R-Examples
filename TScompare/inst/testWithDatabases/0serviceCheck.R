Sys.info()

if(!identical(as.logical(Sys.getenv("_R_CHECK_HAVE_MYSQL_")), TRUE)) {
   cat("MYSQL not available. Skipping tests.\n")
   cat("_R_CHECK_HAVE_MYSQL_ setting ", service1, "\n")
 }else if(!identical(as.logical(Sys.getenv("_R_CHECK_HAVE_SQLITE_")), TRUE)) {
   cat("SQLLITE  not available. Skipping tests.\n")
   cat("_R_CHECK_HAVE_SQLLITE_ setting ", service2, "\n")
 } else  
   {
   require("TScompare")
   require("TSMySQL")
   require("TSSQLite")

   user1    <- Sys.getenv("MYSQL_USER")
   cat("user1 set to:", user1, " by env variable MYSQL_USER\n") 
   if ("" != user1) {
       # specifying host as NULL or "localhost" results in a socket connection
       host1    <- Sys.getenv("MYSQL_HOST")
       if ("" == host1)     host1 <- Sys.info()["nodename"] 
       passwd1  <- Sys.getenv("MYSQL_PASSWD")
       if ("" == passwd1)   passwd1 <- NULL
       }
    
   cat("****************************************************************\n")
   cat("WARNING: THIS OVERWRITES TABLES IN TEST DATABASE ON MySQL SERVER\n")
   cat("****************************************************************\n")

   if ("" != user1)
        setup <- try(RMySQL::dbConnect(RMySQL::MySQL(), dbname="test", 
	                 username=user1, password=passwd1, host=host1)) 
   else setup <- try(RMySQL::dbConnect(RMySQL::MySQL(), dbname="test"))#.my.cnf
 
   if (inherits(setup, "try-error")) stop("dbConnect to MySQL db test failed./n")
 
   TSsql::removeTSdbTables(setup, yesIknowWhatIamDoing=TRUE)
   TSsql::createTSdbTables(setup, index=FALSE)

   DBI::dbListTables(setup) 
 
    # test tables are cleaned up in guideCheck.R
    dbDisconnect(setup)

   if ("" != user1)
         con1 <- try(TSconnect("MySQL", dbname="test", 
	                 username=user1, password=passwd1, host=host1)) 
   else  con1 <- try(TSconnect("MySQL", dbname="test"))#user/passwd/host .my.cnf


   if (inherits(con1, "try-error")) {
        TSsql::removeTSdbTables(con, yesIknowWhatIamDoing=TRUE)
        dbDisconnect(con)
        stop("TSconnect to MySQL db test failed./n")
	}

   cat("*******************************************************************\n")
   cat("* WARNING: THIS OVERWRITES TABLES IN TEST DATABASE ON RSQLite SERVER\n")
   cat("*******************************************************************\n")

   # need to cleanup file "test" created by SQLLite but it is 
   # used in guideCheck.R, so cleanup there.
   setup <- try(RSQLite::dbConnect(RSQLite::SQLite(),
                       dbname="TScompareSQLiteTestDB")) # no user/passwd/host
 
   if (inherits(setup, "try-error"))
           stop("dbConnect to SQLite db test failed./n")

   DBI::dbListTables(setup) 

   TSsql::removeTSdbTables(setup, yesIknowWhatIamDoing=TRUE)
   TSsql::createTSdbTables(setup, index=FALSE)
   
   DBI::dbListTables(setup) 
   dbDisconnect(setup)

   con2 <- try(TSconnect("SQLite", dbname="TScompareSQLiteTestDB")) # no user/passwd/host
   if (inherits(con2, "try-error"))
           stop("TSconnect to SQLite db test failed./n")

#   user2    <- Sys.getenv("POSTGRES_USER")
#   cat("user2 set to:", user2, " by env variable POSTGRES_USER\n") 
#   
#   host2    <- Sys.getenv("POSTGRES_HOST")
#   cat("host2 set to:", host2, " by env variable POSTGRES_HOST\n") 
#     } else  {
#       cat("user2 set to empty string.\n") 
#       #(postgres driver may also use PGDATABASE, PGHOST, PGPORT, PGUSER )
#       # The Postgress documentation seems to suggest that it should be
#       #   possible to get the host from the .pgpass file too, but I cannot.
#       
#       if ("" == host2) {
#          host2 <- Sys.getenv("PGHOST")
#          cat("host2 reset to:", host2, "by env variable PGHOST\n") 
#	  } 
#       if ("" == host2) {
#          host2 <- "localhost"  #Sys.info()["nodename"]
#          cat("host2 reset to:", host2, "\n") 
#	  } 
#       con2 <- try(TSconnect("PostgreSQL", dbname="test"))#user/pwd/host in cfg #
#       if (inherits(con2, "try-error"))
#           stop("connection to PostgreSQL db test failed./n")
#       }

   
   cat("*** TSconnect connections established. Disconnecting\n")
  
   dbDisconnect(con1)
   dbDisconnect(con2)
}

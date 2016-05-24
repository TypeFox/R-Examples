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
   require("TSmisc")  

   user1    <- Sys.getenv("MYSQL_USER")
   cat("user1 set to:", user1, " by env variable MYSQL_USER\n") 
   if ("" != user1) {
       # specifying host as NULL or "localhost" results in a socket connection
       host1    <- Sys.getenv("MYSQL_HOST")
       if ("" == host1)     host1 <- Sys.info()["nodename"] 
       passwd1  <- Sys.getenv("MYSQL_PASSWD")
       if ("" == passwd1)   passwd1 <- NULL
       }
    
   cat("Assuming CLEAN test db set up on MySQL and SQLite by 0serviceCheck.R\n")
 
   if ("" != user1) con1 <- try(TSconnect("MySQL", dbname="test", 
	                 username=user1, password=passwd1, host=host1)) 
   else  con1 <- try(TSconnect("MySQL", dbname="test"))#user/passwd/host .my.cnf
   
   if (inherits(con1, "try-error"))
           stop("TSconnect to MySQL db test failed./n")
 

   con2 <- try(TSconnect("SQLite", dbname="TScompareSQLiteTestDB")) # no user/passwd/host
 
   if (inherits(con2, "try-error"))
           stop("TSconnect to SQLite db TScompareSQLiteTestDB failed./n")
   
   cat("*** TSconnect connections established.\n")
   yahoo <- TSconnect("histQuote", dbname="yahoo") 
   #x <- TSget("^ftse", yahoo)
   # yahoo is slow sometimes
   require("zoo")
   require("tframePlus")
   x <- zoo(rnorm(200), order.by=as.Date("2012-01-01") + 1:200)
   seriesNames(x) <- "ftse"
   TSreplace(x, serIDs="ftse", Table="B", con=con1)
   TSreplace(x, serIDs="ftse", Table="B", con=con2)

   #x <- TSget("^gspc", yahoo)
   x <- zoo(rnorm(200), order.by=as.Date("2012-01-01") + 1:200)
   seriesNames(x) <- "gspc"
   TSreplace(x,  serIDs="gspc", Table="B", con=con1)
   TSreplace(x,  serIDs="gspc", Table="B", con=con2)

   #x <- TSget("ibm", con=yahoo, quote = c("Close", "Vol"))
   x <- zoo(matrix(rnorm(400),200,2), order.by=as.Date("2012-01-01") + 1:200)
   seriesNames(x) <- c("ibmClose", "ibmVol")
   TSreplace(x, serIDs=c("ibmClose", "ibmVol"), Table="B", con=con1)
   TSreplace(x, serIDs=c("ibmC",     "ibmV"),	Table="B", con=con2)

   ids <- AllIds(con1)
   print(ids)
   ids <- cbind(ids, ids)

   eq	<- TScompare(ids, con1, con2, na.rm=FALSE)
   summary(eq)

   eqrm <- TScompare(ids, con1, con2, na.rm=TRUE)
   print(summary(eqrm))

   ids <- matrix(c("ftse","gspc","ibmClose", "ibmVol",
 		   "ftse","gspc","ibmC", "ibmV"),4,2)

   print(ids)
   eq	<- TScompare(ids, con1, con2, na.rm=FALSE)
   print(summary(eq))

   eqrm <- TScompare(ids, con1, con2, na.rm=TRUE)
   print(summary(eqrm))

   cat("*** Remove TSdbi tables in MySQL test db\n")
   require("TSsql")
   removeTSdbTables(con1, yesIknowWhatIamDoing=TRUE)

   cat("*** Disconnecting.\n")
   dbDisconnect(con1)
   dbDisconnect(con2)

   # cleanup file created by SQLLite in 0serviceCheck.R
   unlink("TScompareSQLiteTestDB")

}

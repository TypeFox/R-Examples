Sys.info()

if(!identical(as.logical(Sys.getenv("_R_CHECK_HAVE_MYSQL_")), TRUE)) {
   cat("MYSQL not available. Skipping tests.\n")
   cat("_R_CHECK_HAVE_MYSQL_ setting ", service1, "\n")
 } else {
   require("TScompare")
   require("TSMySQL")

   user1    <- Sys.getenv("MYSQL_USER")
   if ("" != user1) {
       # specifying host as NULL or "localhost" results in a socket connection
       host1    <- Sys.getenv("MYSQL_HOST")
       if ("" == host1)     host1 <- Sys.info()["nodename"] 
       passwd1  <- Sys.getenv("MYSQL_PASSWD")
       if ("" == passwd1)   passwd1 <- NULL
     }

   if ("" != user1){
      con1 <- try(TSconnect("MySQL", dbname="wfsv", 
	        username=user1, password=passwd1, host=host1), silent=TRUE) 
      con2 <- try(TSconnect("MySQL", dbname="ets", 
	        username=user1, password=passwd1, host=host1), silent=TRUE) 
   }else {
      con1 <- try(TSconnect("MySQL", dbname="wfsv"), silent=TRUE)#use~/.my.cnf
      con2 <- try(TSconnect("MySQL", dbname="ets"),  silent=TRUE)#use~/.my.cnf
      }


   if      (inherits(con1, "try-error"))
           cat("wfs connection not available. Skipping 0wfsCheck tests.")
   else if (inherits(con2, "try-error"))
           cat("ets connection not available. Skipping 0wfsCheck tests.")
   else {
      ids <- AllIds(con1)[1:100] #4677 in 2011-05-13
      if(!is.null(AllPanels(con1)))
              stop("Bad result. wfsv does not have panels.")
      if( is.null(AllVintages(con1))) 
              stop("Bad result. wfsv has vintages.")

      ids <- cbind(ids, ids)
      eq   <- TScompare(ids, con1, con2, na.rm=FALSE)
      print(summary(eq))
      eqrm <- TScompare(ids, con1, con2, na.rm=TRUE)
      print(summary(eqrm))
  
      cat("**************        disconnecting.\n")
      dbDisconnect(con1)
      dbDisconnect(con2)
      }
}

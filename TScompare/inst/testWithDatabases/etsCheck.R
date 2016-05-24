if (FALSE){# disable until padi can be replaced by json

service1 <- Sys.getenv("_R_CHECK_HAVE_MYSQL_")
service2 <- Sys.getenv("_R_CHECK_HAVE_PADI_")

# save.image("etsCheck0.RData")

if(!identical(as.logical(service1), TRUE) |
   !identical(as.logical(service2), TRUE)) {
   cat("PADI or MYSQL not available. Skipping tests.\n")
   cat("_R_CHECK_HAVE_MYSQL_ setting ", service1, "\n")
   cat("_R_CHECK_HAVE_PADI_ setting ", service2, "\n")
} else  {

   require("TScompare")
   require("TSMySQL")
   require("TSpadi")

   user1    <- Sys.getenv("MYSQL_USER")
   if ("" != user1) {
       # specifying host as NULL or "localhost" results in a socket connection
       host1    <- Sys.getenv("MYSQL_HOST")
       if ("" == host1)     host1 <- Sys.info()["nodename"] 
       passwd1  <- Sys.getenv("MYSQL_PASSWD")
       if ("" == passwd1)   passwd1 <- NULL
     }

   con1 <- if ("" != user1)  
            tryCatch(TSconnect("MySQL", dbname="ets", username=user1, password=passwd1, host=host1)) 
      else  tryCatch(TSconnect("MySQL", dbname="ets")) # pass user/passwd/host in ~/.my.cnf


   user2    <- Sys.getenv("PADI_USER")
   if ("" != user2) {
       # specifying host as NULL or "localhost" results in a socket connection
       host2    <- Sys.getenv("MYSQL_HOST")
       if ("" == host2)     host2 <- Sys.info()["nodename"] 
       passwd2  <- Sys.getenv("MYSQL_PASSWD")
       if ("" == passwd2)   passwd2 <- NULL
     }
    
   con2 <-  if ("" != user2) 
            try(TSconnect("padi", dbname="ets", username=user, password=passwd, host=host2)) 
       else try(TSconnect("padi", dbname="ets")) # pass user/passwd/host in ~/.padi.cfg

   if (!inherits(con1, "try-error") & !inherits(con2, "try-error")) {
      ids <- AllIds(con1)
      if(!is.null(AllPanels(con1)))   stop("Bad result. ets does not have panels.")
      if(!is.null(AllVintages(con1))) stop("Bad result. ets does not have vintages.")
      ids <- cbind(ids, ids)
      if (FALSE) {
        ids13 <- ids[20001:30000,]
        ids15 <- ids[40001:50000,]

        #save.image("etsCheck0b.RData")
        # This had errors because zoo tf returned by sql causes failure when 
	#  compare with ts tf returned by padi for weekly data: 
	#    Error in if (tf[3] != fr) stop("frequencies must be that same.") :
	#  These are now 
	#  returned as FALSE in tfwindow comparison.
	eq3x   <- TScompare(ids13, con1, con2, na.rm=FALSE) # errors
        print(summary(eq3x))
        eq5   <- TScompare(ids15, con1, con2, na.rm=FALSE) #warnings to check
        print(summary(eq5))

        #save.image("etsCheck0c.RData")
	}

      ids <- ids[1:10000,] # nice to do all, but too slow for regular build
     
      eq   <- TScompare(ids, con1, con2, na.rm=FALSE)
      print(summary(eq))
      eqrm <- TScompare(ids, con1, con2, na.rm=TRUE)
      print(summary(eqrm))
  
      cat("**************        disconnecting ets\n")
      dbDisconnect(con1)
      #dbDisconnect(con2)
      }
}
}

# test that the fame package is working with a server

#  (at BoC I need environment  export FAME=/apps/fame92r2 ) for fame
# export FAME=/apps/fame92r2

if(identical(as.logical(Sys.getenv("_R_CHECK_HAVE_FAME_")), TRUE)) {

 # require("tis")
 require("fame")

cat("***********   test fame server (ets)\n")
# this is using Fame drivers to talk to local or remote dbs.

  if(!fameRunning()) fameStart(workingDB = FALSE)
 
  Id <- try(fameDbOpen("ets /home/ets/db/etsmfacansim.db", accessMode = "read"))
  if(inherits(Id, "try-error") ) 
      stop("Could not establish fameConnection to ets /home/ets/db/etsmfacansim.db")

  fameDbClose(Id) #this needs fix in fame 2.11

cat("***********   reading from ets /home/ets/db/etsmfacansim.db\n")

   r <- getfame("V122646", db="ets /home/ets/db/etsmfacansim.db", save = FALSE, 
             envir = parent.frame(),
             start = NULL, end = NULL, getDoc = FALSE)[[1]]

cat("***********   reading from ets /home/ets/db/etsmfacansim.db using con\n")
# this is using  connection which uses remote Fame server.
 
  # requires fame 2.11 
  #     user = "", password = ""
  con <- fameConnection(service = "2959", host = "ets", stopOnFail = TRUE)

   r2 <- getfame("V122646", db="/home/ets/db/etsmfacansim.db",
      connection=con, save = FALSE, 
             envir = parent.frame(),
             start = NULL, end = NULL, getDoc = FALSE)[[1]]
 
# check that connection was not closed
   r3 <- getfame("V122646", db="/home/ets/db/etsmfacansim.db",
      connection=con, save = FALSE, 
             envir = parent.frame(),
             start = NULL, end = NULL, getDoc = FALSE)[[1]]
 
 if (1e-14 < max(abs(r - r2))) stop("Fame call results differ.") 

} else  cat("FAME not available. Skipping tests.")

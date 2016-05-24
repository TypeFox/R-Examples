# test that the fame package is working with Fame Server

#  (at BoC I need environment  export FAME=/apps/fame92r2 ) for fame
# export FAME=/apps/fame92r2

if(identical(as.logical(Sys.getenv("_R_CHECK_HAVE_FAME_")), TRUE)) {

cat("**************   test fame Server\n")

  require("TSfame")

cat("**************        connecting ets fame server\n")
# this is using (fake) connection which uses Fame drivers to talk to
# local or remote dbs.

con <- TSconnect("fame", dbname="ets /home/ets/db/etsmfacansim.db", "read") 

if(!inherits(con, "try-error") ) {
   z <- TSget("V122646", con=con)

   if(any(start(z) != c(1969,1))) stop("Error reading V122646.")
   tfplot(z)

   dbDisconnect(con)
   } else  cat("ets fame server not available. Skipping tests.")

cat("***********   reading from ets /home/ets/db/etsmfacansim.db using con\n")
# this is using  connection which uses remote Fame server.

con2 <- TSconnect("fameServer", dbname="/home/ets/db/etsmfacansim.db", 
    service = "2959", host = "ets", stopOnFail = TRUE) 

r <- TSget("V122646", con=con2)
r <- TSget(c("V122646", "B171"), con=con2)

TSdates("V122646", con=con2)
TSdoc("V122646", con=con2) 
TSdescription("V122646", con=con2) 
TSlabel("V122646", con=con2) 

dbDisconnect(con2)

} else  cat("FAME not available. Skipping tests.")

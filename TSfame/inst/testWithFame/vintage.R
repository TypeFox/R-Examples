
if(identical(as.logical(Sys.getenv("_R_CHECK_HAVE_FAME_")), TRUE)) {

require("TSfame")
require("tfplot")

cat("**************  vintage  examples with fame Server\n")

cat("**************   connecting ets fame server using local Fame server\n")
# this is using (fake) connection which uses Fame drivers to talk to
# local or remote dbs.

# alternatively the vintage names (dates) and files paths can be in a file, and
#  us vintageMap() to read it
dbs <- paste("ets /home/ets5/mfadata/etsmfacansim_", c(
             "20110513.db", "20060526.db", "20110520.db"), sep="")
names(dbs) <- c("2011-05-13", "2006-05-26", "2011-05-20")
	     
con <- TSconnect("fame", dbname=dbs, "read") 

if(!inherits(con, "try-error") ) {   
   z <- TSget("V122646", con=con, vintage=c("2011-05-13", "2006-05-26"))

   if(any(start(z) != c(1969,1))) stop("Error reading V122646.")
   tfplot(z)

   dbDisconnect(con)
   }
	     
con <- TSconnect("fame", dbname=dbs, "read", current="2011-05-13") 

if(!inherits(con, "try-error") ) {
   z <- TSget("V122646", con=con) 
   z <- TSget(c("V122646", "V122647"),con=con) 
   
   if(any(start(z) != c(1969,1))) stop("Error reading V122646.")
   tfplot(z)

   dbDisconnect(con)
   }

cat("***********   connecting ets fame server using remote Fame server\n")
# this is using  connection which uses remote Fame server.

dbs <- paste("/home/ets5/mfadata/etsmfacansim_", c(
             "20110513.db", "20060526.db", "20110520.db"), sep="")
names(dbs) <- c("2011-05-13", "2006-05-26", "2011-05-20")

con2 <- TSconnect("fameServer", dbname=dbs, 
    service = "2959", host = "ets", stopOnFail = TRUE, current="2011-05-13") 

if(!inherits(con2, "try-error") ) {
   r <- TSget("V122646", con=con2, vintage=c("2011-05-13", "2006-05-26"))
   TSdates("V122646", con=con2, vintage=c("2006-05-26"))

   # next just gets current by default
   r <- TSget(c("V122646", "V122647"), con=con2)

   TSdates("V122646", con=con2)
   TSdoc("V122646", con=con2) 
   TSdescription("V122646", con=con2) 
   TSlabel("V122646", con=con2) 

   dbDisconnect(con2)
   }
} else  cat("FAME not available. Skipping tests.")

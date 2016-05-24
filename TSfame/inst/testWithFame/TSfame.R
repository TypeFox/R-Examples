#  (at BoC I need environment  export FAME=/apps/fame92r2 ) for fame
# export R_LIBS=/home/mfa/gilp/serveR.myBase/gcc-3.4.6/R-2.7.1 
# export FAME=/apps/fame92r2

require("tfplot")

if(identical(as.logical(Sys.getenv("_R_CHECK_HAVE_FAME_")), TRUE)) {

require("TSfame")

cat("**************        connecting testFame.db\n")

con <- TSconnect("fame", dbname="testFame.db") 

cat("**************        examples\n")

seriesA  <- tis(1:24,	   start = c(2002, 1), freq = 12)
seriesNames(seriesA) <- "seriesA"
seriesA2 <- tis(rnorm(24), start = c(2002, 1), freq = 12)
seriesNames(seriesA2) <- "seriesA2"

seriesB <- tis(rnorm(104), start = c(2002, 1), tif = "wmonday")
seriesNames(seriesB) <- "seriesB"

seriesC <- tis(rnorm(104), start = c(2002, 1), tif = "wfriday")
seriesNames(seriesC) <- "seriesC"
start(seriesC)


TSdoc(seriesB) <- paste("Line", 1:4, "of seriesB documentation", collapse="\n")
TSdescription(seriesB) <- "seriesB description"

cat("*******  cleanup to db\n")

if( TSexists("seriesA", con)) TSdelete("seriesA", con)
if( TSexists("seriesA2", con)) TSdelete("seriesA2", con)
if( TSexists("seriesB", con)) TSdelete("seriesB", con)
if( TSexists("seriesC", con)) TSdelete("seriesC", con)
if( TSexists("tsseries", con)) TSdelete("tsseries", con)
if( TSexists("tsseriesM1", con)) TSdelete("tsseriesM1", con)
if( TSexists("tsseriesM2", con)) TSdelete("tsseriesM2", con)

if( TSexists(c("tsseriesM1", "tsseriesM2"), con)) 
    TSdelete(c("tsseriesM1", "tsseriesM2"), con)

cat("*******  write to db\n")

if(!TSput(seriesA,  con))    stop("TSput seriesA error.")
if(!TSput(seriesA2, con))    stop("TSput seriesA2 error.")
if(!TSput(seriesB, con))     stop("TSput seriesB error.")
if(!TSreplace(seriesB, con)) stop("TSreplace seriesB error.")

tsseries <- ts(rnorm(30), start=c(1999,2), frequency=12)
seriesNames(tsseries) <- "tsseries"
if(!TSput(tsseries,  con))    stop("TSput tsseries error.")

tsseriesM <- ts(matrix(rnorm(60), 30,2), start=c(1999,2), frequency=12)
seriesNames(tsseriesM) <- c("tsseriesM1","tsseriesM2")
if(!TSput(tsseriesM,  con))    stop("TSput tsseriesM error.")

cat("*******  read from db\n")

z  <- TSget("tsseries", con)
if(max(abs(c(z) - c(tsseries))) > 1e-15) stop("TSget tsseries error.")

if(any(start(z) != c(1999,2))) stop("TSget tsseries start error.")

if(any(start(TSdates("tsseries", con))[[1]] != c(1999,2)))
     stop("TSdates tsseries start error.")


z  <- TSget(c("tsseriesM1","tsseriesM2"), con)
if(max(abs(c(z) - c(tsseriesM))) > 1e-15) stop("TSget tsseriesM error.")

if(any(start(z) != c(1999,2))) stop("TSget tsseriesM start error.")

if(any(start(TSdates(c("tsseriesM1","tsseriesM2"), con))[[1]] != c(1999,2)))
if(any(start(TSdates(seriesNames(tsseriesM), con))[[1]] != c(1999,2)))
     stop("TSdates tsseriesM start error.")


a  <- TSget("seriesA", con)
if(max(abs(c(a) - c(seriesA))) > 1e-15) stop("TSget seriesA error.")

if(any(start(TSdates("seriesA", con))[[1]] != c(2002, 1)))
     stop("TSget seriesA start error.")

if(!TSexists("seriesA", con)) stop("TSget seriesA exists error.")
if( TSexists("seriesZ", con)) stop("TSget seriesZ !exists error.")

a2 <- TSget("seriesA2", con)
if(max(abs(c(a2) - c(seriesA2))) > 1e-15) stop("TSget seriesA2 error.")


z <- TSget(c("seriesA","seriesA2"), con)
tfplot(z)

b  <- TSget("seriesB", con, TSrepresentation="zoo")
if(max(abs(c(b) - c(seriesB))) > 1e-15) stop("TSget seriesB error.")

if(as.Date(start(b))!= "2002-01-07")              stop("seriesB start date error.")
if(weekdays(start(seriesB)) != weekdays(start(b)))stop("seriesB weekdays error.")
if(weekdays(start(b)) != "Monday")                stop("seriesB start weekdays error.")
if(weekdays(start(b)) != weekdays(end(b)))        stop("seriesB weekdays error.")

tfplot(b)

tframe(seriesB)

tfwindow(seriesB, tf=NULL)
bi <- TSget("seriesB", con, TSrepresentation="tis")
if(max(abs(c(bi) - c(seriesB))) > 1e-15) stop("TSget bi error.")

if(base::as.Date(start(bi))!= "2002-01-07")        stop("bi start date error.")
if(weekdays(start(seriesB)) != weekdays(start(bi)))stop("bi weekdays error.")
if(weekdays(start(bi)) != "Monday")                stop("bi start weekdays error.")
if(weekdays(start(bi)) != weekdays(end(bi)))       stop("bi weekdays error.")

tfplot(bi)

cat("*******  documentation from db\n")

if(TSdoc("seriesB", con) != TSdoc(seriesB)) stop("TSdoc error.")

if(TSdescription("seriesB", con) != TSdescription(seriesB))
           stop("TSdescription error.")

cat("*******  TSdates from db\n")

TSdates("seriesB", con)
TSdates(c("seriesA","seriesA2","seriesB"), con)

cat("*******  delete from db\n")

TSdelete("seriesB", con)

if(!TSput(seriesC, con))     stop("TSput seriesC error.")
if( TSput(seriesC, con))     stop("TSput seriesC error overwriting.")

ci  <- TSget("seriesC", con, TSrepresentation="tis")
start(ci)
start(bi)

cat("*******  timeSeries representation\n")
ci  <- TSget("seriesC", con, TSrepresentation="timeSeries")
if("timeSeries" != class(ci)) stop("timeSeries class object not returned.")
start(ci)
start(bi)

cat("*******  misc\n")

TSdates("seriesA", con)
if(base::as.Date(start(seriesC))!= "2002-01-04") stop("seriesC start date error.")

if(weekdays(start(seriesC)) != "Friday") stop("seriesC start weekdays error.")
if(weekdays(start(seriesC)) != weekdays(end(seriesC))) stop("seriesC weekdays error.")


#dayOfWeek(start(seriesC))

cat("**************        disconnecting test\n")
dbDisconnect(con)

} else  cat("FAME not available. Skipping tests.")

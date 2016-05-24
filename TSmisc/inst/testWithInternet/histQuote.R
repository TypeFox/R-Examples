require("TSmisc")
require("tfplot")
cat("************** TShistQuote  Examples ******************************\n")

con <- TSconnect("histQuote", dbname="yahoo") 

x <- TSget("^ftse", con)  #"^gdax" previously but causing problems
plot(x)
tfplot(x)

# debugging info
str(x)
TSrefperiod(x)

if ("Close" != TSrefperiod(x)) stop("TSrefperiod error, test 1.") 
TSdescription(x) 

# Next produces warnings when yahoo messes up dates and puts the previous 
# previous day's dates on holidays: "some methods for zoo objects do not 
# work if the index entries in order.by are not unique" 
# The warnings can be suppressed
#  x2 <- suppressWarnings( TSget("^gspc", con))
# but an error will occur further below when two series are retrieved and
# zoo is used to merge them.

x2 <- TSget("^gspc", con)
tfplot(x2)
plot(x2)
if ("Close" != TSrefperiod(x2)) stop("TSrefperiod error, test 2.") 
TSdescription(x2) 

options(TSconnection=con)

x <- TSget(c("^ftse","^gspc"), con=con)
plot(x)
tfplot(x)
if (! all(TSrefperiod(x) %in% "Close")) stop("TSrefperiod error, test 3.")
TSdescription(x) 

x <- TSget("ibm", quote = c("Close", "Vol"), TSrepresentation="timeSeries")
if("timeSeries" != class(x)) stop("timeSeries class object not returned.")
plot(x)
tfplot(x)
if(!all(TSrefperiod(x) == c("Close", "Vol"))) stop("TSrefperiod error, test 4.")
TSdescription(x) 

tfplot(x, xlab = TSdescription(x))
tfplot(x, Title="IBM", start="2007-01-01")

conO <- TSconnect("histQuote", dbname="oanda") 

# oanda has max 500 data points
z <- TSget("EUR/USD", conO, start=Sys.Date() - 495)
plot(z)
tfplot(z)
if ("Close" != TSrefperiod(z)) stop("TSrefperiod error, test 5.") 

tfplot(z, Title = "EUR/USD")
tfplot(z, Title = "EUR/USD", start="2007-01-01")
tfplot(z, Title = "EUR/USD", start="2007-03-01")
tfplot(z, Title = "EUR/USD", start=Sys.Date()-14, end=Sys.Date(),
   xlab = format(Sys.Date(), "%Y"))


TSdates(c("^ftse","^gspc", "ibm"), con) # note default start

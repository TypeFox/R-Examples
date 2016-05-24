require("TSsdmx")
require("tframe")

# RJSDMX::sdmxHelp()  # can be useful for finding series identifiers, etc

eurostat <- TSconnect("sdmx", dbname="EUROSTAT")

# this worked in Dec 2014 but the series seems to have disappeared
# z <- TSget("ei_nama_q.Q.MIO-EUR.NSA.CLV2000.NA-B1G.IT", eurostat)

z <- TSget("ei_nama_q.Q.MIO-EUR.SWDA.CP.NA-P72.IT", eurostat)

if (seriesNames(z) != "ei_nama_q.Q.MIO-EUR.SWDA.CP.NA-P72.IT")
    stop("seriesNames not set properly in eurostat test 1.")
    
TSmeta(z)

if (! all(c(1980, 1) == start(z))) stop("eurostat test 1 start date has changed.")


z <- TSget("ei_nama_q.Q.MIO-EUR.SWDA.CP.NA-P72.IT",
           start="1990-Q1", end="2012-Q2", eurostat)

if (! all(c(1990, 1) == start(z))) stop("eurostat test 2 start date failure.")
if (! all(c(2012, 2) ==   end(z))) stop("eurostat test 2  end  date failure.")


#z <- TSget('ei_nama_q.Q.MIO-EUR.NSA.CP.*.IT', eurostat) # all NaN

# at one time this had 28 series, 23 with data
#z <-  TSget("ei_nama_q.Q.MIO-EUR.NSA.CLV2000.*.IT", eurostat) 
z <-  TSget("ei_nama_q.Q.MIO-EUR.NSA.CP.*.IT", eurostat) 

if (28 != sum(hasData(z, quiet=TRUE)))    stop("eurostat hasData test 1 changed.") # previously 23
if (28 != length(hasData(z, quiet=TRUE))) stop("eurostat hasData test 2 changed.")

hasDataCount(z)
hasDataNames(z)

hasDataDescriptions(z)

#  This is a useful check to know if a series has data
if (! ("ei_nama_q.Q.MIO-EUR.NSA.CP.NA-P72.IT"
                 %in% hasDataNames(z))) stop("eurostat hasData test 3 changed.")

##  vector of serIDs

z <-  TSget(c("ei_nama_q.Q.MIO-EUR.NSA.CP.NA-P7.IT",        
              "ei_nama_q.Q.MIO-EUR.NSA.CP.NA-P71.IT",        
              "ei_nama_q.Q.MIO-EUR.NSA.CP.NA-P72.IT"),
	    start="1990-Q1", end="2012-Q2",eurostat) 

if (! all(c(1990, 1) == start(z))) stop("eurostat vector test 1 start date failure.")
if (! all(c(2012, 2) ==   end(z))) stop("eurostat vector test 1  end  date failure.")
if ( 4 !=  frequency(z)) stop("eurostat vector test 1  frequency  date failure.")


z <-  TSget(c("ei_nama_q.Q.MIO-EUR.NSA.CP.NA-P7.IT",	    
              "ei_nama_q.Q.MIO-EUR.NSA.CP.NA-P71.IT",	     
              "ei_nama_q.Q.MIO-EUR.NSA.CP.NA-P72.IT"),
	    start=c(1990,1), end=c(2012,2), eurostat) 

if (! all(c(1990, 1) == start(z))) stop("eurostat vector test 1 start date failure.")
if (! all(c(2012, 2) ==   end(z))) stop("eurostat vector test 1  end  date failure.")
if ( 4 !=  frequency(z)) stop("eurostat vector test 1  frequency  date failure.")

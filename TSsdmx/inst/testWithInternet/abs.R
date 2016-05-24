########################## Australian Bureau of Statistics ####################
# http://www.abs.gov.au

require("TSsdmx")
require("tframe")

abs <- TSconnect("sdmx", dbname="ABS")

##################  quarterly #################

z <- TSget("BOP.1.100.10.Q", start=c(1990, 1), end=c(2012, 2), abs) 

if (seriesNames(z) != "BOP.1.100.10.Q")
    stop("seriesNames not set properly in abs test 1.")

if (! all(c(1990, 1) == start(z))) stop("abs test 1 start date failure.")
if (! all(c(2012, 2) ==   end(z))) stop("abs test 1  end  date failure.")
if   (4 != frequency(z))           stop("abs test 1 frequency failure.")

##################  monthly #################

z <- TSget('RT.0.2.15.10.M', start=c(1990, 1), end=c(2012, 2), names="Retail Trade", abs) 

if ("Retail Trade" != seriesNames(z))    stop("abs test 2 seriesName failure.")
if (! all(c(1990, 1) == start(z))) stop("abs test 2 start date failure.")
if (! all(c(2012, 2) == end(z)))   stop("abs test 2  end  date failure.")
if   (12 != frequency(z))          stop("abs test 2 frequency failure.")

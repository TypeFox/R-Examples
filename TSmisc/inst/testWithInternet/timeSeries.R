require("TSmisc")
require("timeSeries")
require("tfplot")

##################################################
##################################################

#### Data from PiTrading  ########
## http://pitrading.com/free_market_data.htm # free futures data 
## http://pitrading.com/free_eod_data/INDU.zip
##################################################
##################################################
  pit <- TSconnect("zip", dbname="http://pitrading.com/free_eod_data")


  z <- TSget("INDU", pit, TSrepresentation="timeSeries")
  tfplot(z)

  options(TSrepresentation="timeSeries")

  z <- TSget(c("EURUSD", "GBPUSD"), pit)
  if("timeSeries" != class(z)) stop("timeSeries class object not returned.")
  tfplot(z)
 
  TSrefperiod(z) 
  TSdescription(z) 

  z <- TSget(c("AD", "CD"), pit, quote="Close")

  zz <- window(z, start=as.timeDate("2007-01-01"), end=end(z))
  zz <- window(z, start="2007-01-01", end=end(z))
  zz <- tfwindow(z, start="2007-01-01", end=end(z))
  zz <- tfwindow(z, start="2007-01-01")
  

  tfplot(z, start="2007-01-01",
         Title="Australian and Canadian Dollar Continuous Contract, Close")

  unlink("Rplots.pdf")

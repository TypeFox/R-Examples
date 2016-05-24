require("TSmisc")
require("xts")
require("tfplot")

##################################################
##################################################

#### Data from PiTrading  ########
## http://pitrading.com/free_market_data.htm # free futures data 
## http://pitrading.com/free_eod_data/INDU.zip
##################################################
##################################################
  pit <- TSconnect("zip", dbname="http://pitrading.com/free_eod_data")

  z <- TSget("INDU", pit, TSrepresentation=xts)
  tfplot(z)

  options(TSrepresentation=xts)

  z <- TSget(c("EURUSD", "GBPUSD"), pit)
  tfplot(z)
 
  TSrefperiod(z) 
  TSdescription(z) 

  z <- TSget(c("AD", "CD"), pit, quote="Close")
  tfplot(z, start="2007-01-01",
         Title="Australian and Canadian Dollar Continuous Contract, Close")

  unlink("Rplots.pdf")

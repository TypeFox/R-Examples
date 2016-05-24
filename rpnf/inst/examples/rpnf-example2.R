#' This example shows how to get a recent plot for any given yahoo symbol
#' 
### Initialize libraries
library(rpnf) # Load rpnf library
library(quantmod) # Load quantmod library

### Define wrapper for quantmod download function
downloadData <- function(symbol="GOOG") {
  stockData <- new.env() #Make a new environment for quantmod to store data in
  startDate = as.Date(Sys.Date()-3*365) #Specify period of time we are interested in
  endDate = as.Date(Sys.Date()) # End date is today
  getSymbols(symbol, env = stockData, src = "yahoo", from = startDate, to = endDate)   #Download the stock history (for all tickers)
  ### extract stock quotes from enviorment, rename them and store to data.frame 
  symbolData <- eval(parse(text=paste("OHLC(stockData$",sub("^\\^","",symbol),")",sep="")))
  names(symbolData)<-c("open","high","low","close")
  symbolData
}

### The example code starts
# Define (yahoo) symbol to be processed
symbol <- "^GDAXI" # e.g. "GOOG", "^DJI", "^GDAXI", "A"
# download stock quotes for last 3 years
data <- downloadData(symbol)
# Generate point and figure informations
log <- TRUE
boxsize <- getLogBoxsize(percent=3)
pnf <- pnfprocessor(high=data$high,low=data$low,date=index(data),boxsize=boxsize,log=log)
# View result data
View(pnf)
# Plot result as text to console
pnfplottxt(pnf,boxsize=boxsize,log=log,main=paste("P&F Plot ",symbol))
# Plot result as graphic
pnfplot(pnf,boxsize=boxsize,log=log,main=paste("P&F Plot ",symbol))

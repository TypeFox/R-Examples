# Script for parsing gold price data series from the gold.org website.

# Peter Carl

# DETAILS
# Parse monthly and daily close prices from the spreadsheet:
# http://www.gold.org/download/value/stats/statistics/xls/gold_prices.xls

# Load needed packages:
require(zoo)
require(gdata)
require(FinancialInstrument)
require(quantmod)
# Set the working directory, where there's a .incoming folder that contains
# the downloaded spreadsheet.
filesroot = "~/Data/gold"

# Create and set the working directory if it doesn't exist
if (!file.exists(filesroot))
  dir.create(filesroot, mode="0777")

# Create and set the .incoming directory if it doesn't exist
if (!file.exists(paste(filesroot, "/.incoming", sep="")))
  dir.create(paste(filesroot, "/.incoming", sep=""), mode="0777")
setwd(paste(filesroot, "/.incoming", sep=""))

# Remove the old file from .incoming
if(file.exists("gold_prices.xls"))
  system(paste("rm gold_prices.xls"))

# Download the xls workbook directly from the web site:
system("wget http://www.gold.org/download/get/value/stats/statistics/xls/gold_prices.xls")

if(!file.exists("gold_prices.xls"))
  stop(paste("No spreadsheet exists.  Download the spreadsheet to be processed from www.gold.org/download/value/stats/statistics/xls into ", filesroot, "/.incoming", sep=""))
  
x = read.xls("gold_prices.xls", sheet="Monthly_EndofPeriod", pattern="USD")
# This can be repeated with sheet="Daily_EndofPeriod" and processed similarly

# Remove the first three columns of NAs
x=x[,-3:-1]
# Remove repeated data
x=x[,c(-10,-19,-21,-23,-24)]

# Get the descriptions to add as attributes
x.attr = x[1,]
x=x[-1,]


# Get attributes and labels
categoryNames = x.attr
symbolNames=paste("GOLD", colnames(x), ".M",sep="")
ISOdates = as.Date(x[,1])

for(i in 2:length(symbolNames)) {
  # check to make sure directories exist for each symbol
  dir.create(paste(filesroot, symbolNames[i], sep="/"), showWarnings = FALSE, 
  recursive = FALSE, mode = "0777")
}

# Parse the columns into individual price objects
for( i in 2:dim(x[,-1])[2]){
  x.xts = as.xts(as.numeric(gsub(",","",x[,2])), order.by=ISOdates)
  R.xts = Return.calculate(x.xts)
  x.xts = cbind(x.xts, R.xts)
  colnames(x.xts)=c("Close", "Returns")
  xtsAttributes(x.xts) <- list(Description = paste(categoryNames[,i], "Prices"))
#   assign(symbolNames[i], x.xts)
  save(x.xts, file=paste(filesroot, symbolNames[i], paste(symbolNames[i], ".rda", sep=""), sep="/"))
  # Describe the metadata for each index
  instrument(symbolNames[i], currency="USD", multiplier=1, tick_size=.01, start_date=head(index(x.xts),1), description=paste(categoryNames[,i], "per troy ounce"), data="CR", source="gold.org", assign_i=TRUE)
}

# Now, whenever you log in you need to register the instruments.  This
# might be a line you put into .Rprofile so that it happens automatically:
# require(quantmod) # this requires a development build after revision 560 or so.
setSymbolLookup.FI(base_dir=filesroot, split_method='common')

# Now you should be able to:
getSymbols("GOLDUSD.M")
chartSeries(Cl(GOLDUSD.M), theme="white")
charts.PerformanceSummary(GOLDUSD.M["2004::","Returns"], ylog=TRUE, wealth.index=TRUE, main = "Gold Monthly Returns")
tail(GOLDUSD.M)
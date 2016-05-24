# Script for downloading and parsing a monthly total return series from
# http://www.mscibarra.com/
#

# Peter Carl

# Currently, this script assumes that a single index is in the spreadsheet.
# For cases where there are multiple indexes, this script will need to be
# extended.

# TODO: Modify to take in four indexes in a particular order

# Load needed packages:
require(xts)
require(PerformanceAnalytics)
require(gdata)
require(FinancialInstrument)
require(quantmod)
currency("USD")

# Set the working directory, where there's a .incoming folder that contains
# the downloaded spreadsheet.
filesroot = "~/Data/MSCI"

# Create and set the working directory if it doesn't exist
if (!file.exists(filesroot))
  dir.create(filesroot, mode="0777")

# Create and set the .incoming and symbol directories if they don't exist
if (!file.exists(paste(filesroot, "/.incoming", sep="")))
  dir.create(paste(filesroot, "/.incoming", sep=""), mode="0777")
setwd(paste(filesroot, "/.incoming", sep=""))

if(!file.exists("historyIndex.xls"))
  stop(paste("No spreadsheet exists.  Download the spreadsheet to be processed from http://www.msci.com into ", filesroot, "/.incoming", sep=""))

# Read the first sheet in the xls workbook directly from the working directory:
x = read.xls("historyIndex.xls", pattern="Date")
x = x[-((dim(x)[1]-15):dim(x)[1]),] # trim off last 16 lines of disclaimer
x.dates = paste(substring(x[,1],1,6),substring(x[,1],7,10)) # unmangle the dates
ISOdates = as.Date(x.dates, format="%b %d %Y")

symbolNames = colnames(x[,-1])
description = paste("MSCI", gsub(pattern = "\\.", replacement = " ", symbolNames))
description = gsub(pattern = "  ", replacement = " ", description)
symbolNames = gsub(pattern = ".Standard.*", replacement = "", symbolNames)
symbolNames = paste("MSCI.", symbolNames, ".M.IDX", sep="")
symbolNames = gsub(pattern = "\\.\\.", replacement = "\\.", symbolNames) # Remove double dots

for(i in 1:length(symbolNames)) {
  # check to make sure directories exist for each symbol
  dir.create(paste(filesroot, symbolNames[i], sep="/"), showWarnings = FALSE, 
             recursive = FALSE, mode = "0777")
}

# Parse the columns into individual price objects
print("Processing columns as symbols...")
for( i in 1:length(symbolNames)){
  # x.prices = as.numeric((sub(",","", x[,-1], fixed=TRUE)))
  x.xts = as.xts(as.numeric(sub(",","", x[,i+1], fixed=TRUE)), order.by=ISOdates)
  R.xts = Return.calculate(x.xts)
  x.xts = cbind(x.xts, R.xts)
  index(x.xts) = as.Date(as.yearmon(index(x.xts)), frac=1) # resets date to last day of the month rather than last business day
  colnames(x.xts)=c(paste(symbolNames[i],".Close", sep=""), paste(symbolNames[i], ".Returns", sep=""))
  xtsAttributes(x.xts) <- list(Description = paste(description[i], "Index"))
  save(x.xts, file=paste(filesroot, symbolNames[i], paste(symbolNames[i], ".rda", sep=""), sep="/"))
  # Describe the metadata for each index
  instrument(symbolNames[i], currency="USD", multiplier=1, tick_size=.01, start_date=head(index(x.xts),1), description=paste(description[i], "Index", sep=""), data="CR", source="MSCI", frequency="Monthly", assign_i=TRUE)
  print(paste(symbolNames[i], description[i]))
}

# Now, whenever you log in you need to register the instruments.  This
# might be a line you put into .Rprofile so that it happens automatically:
# require(quantmod) # this requires a development build after revision 560 or so.
setSymbolLookup.FI(base_dir=filesroot, split_method='common')

# Now you should be able to:
getSymbols("MSCI.WORLD.M.IDX")
getSymbols(symbolNames)
chartSeries(Cl(MSCI.WORLD.M.IDX), theme="white")
head(MSCI.WORLD.M.IDX)
# Script for downloading and parsing a monthly total return series of the FTSE NAREIT U.S. Real Estate Index
#   - Download the spreadsheet at:
#     http://returns.reit.com/returns/MonthlyHistoricalReturns.xls
#

# Peter Carl

# Load needed packages:
require(xts)
require(gdata)
require(FinancialInstrument)
require(quantmod)

# Set the working directory, where there's a .incoming folder that contains
# the downloaded spreadsheet.
filesroot = "~/Data/NAREIT"

# Create and set the working directory if it doesn't exist
if (!file.exists(filesroot))
  dir.create(filesroot, mode="0777")
  
if (!file.exists(paste(filesroot, "/NAREIT.M.IDX", sep="")))
  dir.create(paste(filesroot, "/NAREIT.M.IDX", sep=""), mode="0777")
  
# Download the first sheet in the xls workbook directly from the web site:
x = read.xls("http://returns.reit.com/returns/MonthlyHistoricalReturns.xls", pattern="Date", sheet="Index Data", stringsAsFactors=FALSE)

# We'll focus on the first three columns, the date, the total return, and the total return index.  The next columns divide it into price, income and dividend yield components, and after that are the individual sectors.
# @TODO: create symbols for the rest of the columns in the spreadsheet
x.dates = as.Date(as.yearmon(x[,1], format="%b-%y"), frac=1)
x.returns = xts(x[,2]/100, order.by = x.dates)
x.price = as.numeric((sub(",","", x[,3], fixed=TRUE))) # get rid of commas
x.price = xts(x.price, order.by = x.dates)

x.xts = cbind(x.price, x.returns)
colnames(x.xts) = c("NAREIT.M.IDX.Close", "NAREIT.M.IDX.Returns")

# Save it into an rda file on the filesystem
save(x.xts, file=paste(filesroot,"NAREIT.M.IDX/NAREIT.M.IDX.rda", sep="/"))

# Create currencies first:
require(FinancialInstrument)
currency("USD")

# Describe the metadata for the index
instrument("NAREIT.M.IDX", currency="USD", multiplier=1, tick_size=.01, start_date="1971-12-31", description="FTSE NAREIT U.S. Real Estate Index", data="CR", source="reit.com", assign_i=TRUE)

# Now, whenever you log in you need to register the instruments.  This
# might be a line you put into .Rprofile so that it happens automatically:
# require(quantmod) # this requires a development build after revision 560 or so.
setSymbolLookup.FI(base_dir=filesroot, split_method='common')

# Now you should be able to:
getSymbols("NAREIT.M.IDX")
chartSeries(Cl(NAREIT.M.IDX), theme="white")
head(NAREIT.M.IDX)

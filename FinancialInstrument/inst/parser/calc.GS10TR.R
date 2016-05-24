# Use yields of the Constant Maturity 10 year bond series from FRED 
# to calculate total returns.

# Peter Carl

# Originally described by Kenton Russell at his blog TimelyPortfolio,
# posted April 15, 2011
# http://timelyportfolio.blogspot.com/2011/04/historical-sources-of-bond-returns_17.html
#

require(quantmod)
require(PerformanceAnalytics)
require(RQuantLib)

# Set the working directory.
filesroot = "~/Data/FRED"

# Create and set the working directory if it doesn't exist
if (!file.exists(filesroot))
  dir.create(filesroot, mode="0777")
  
if (!file.exists(paste(filesroot, "/GS10TR.IDX", sep="")))
  dir.create(paste(filesroot, "/GS10TR.IDX", sep=""), mode="0777")
  
getSymbols("GS10", src="FRED") #load US Treasury 10y yields from FRED

# Dates should be end of month, not beginning of the month as reported
index(GS10) = as.Date(as.yearmon(index(GS10)), frac=1)

# @TODO: Do this calculation with a longer list of symbols

x.pr <- GS10  #set this up to hold price returns
x.pr[1,1] <- 0
colnames(x.pr) <- "Price Return"
for (i in 1:(NROW(GS10)-1)) {
  x.pr[i+1,1] <- FixedRateBondPriceByYield(yield=GS10[i+1,1]/100, issueDate=Sys.Date(), maturityDate=advance("UnitedStates/GovernmentBond", Sys.Date(), 10, 3), rates=GS10[i,1]/100,period=2)[1]/100-1
}
#total return will be the price return + yield/12 for one month
x.tr <- x.pr + lag(GS10,k=1)/12/100
colnames(x.tr)<-"Total Return"

# Add an index column labeled "Close"
x.idx = 100*cumprod(1+na.omit(x.tr))
x.xts = cbind(x.idx, x.tr)
x.xts[1,1]=100 # base the index at 100
colnames(x.xts) = c("Close", "Return")

# Save it into an rda file on the filesystem
save(x.xts, file=paste(filesroot,"GS10TR.IDX/GS10TR.IDX.rda", sep="/"))

# Create currencies first:
require(FinancialInstrument)
currency("USD")

# Describe the metadata for the index
instrument("GS10TR.IDX", currency="USD", multiplier=1, tick_size=.01, start_date="1953-04-01", description="US 10Y Constant Maturity Total Returns", data="CR", source="fred", assign_i=TRUE)

# Now, whenever you log in you need to register the instruments.  This
# might be a line you put into .Rprofile so that it happens automatically:
# require(quantmod) # this requires a development build after revision 560 or so.
setSymbolLookup.FI(base_dir=filesroot, split_method='common')

# Now you should be able to:
getSymbols("GS10TR.IDX")
chartSeries(Cl(GS10TR.IDX), theme="white")
head(GS10TR.IDX)

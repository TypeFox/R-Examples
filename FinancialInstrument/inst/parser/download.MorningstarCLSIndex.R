# Script for downloading and parsing Morningstar CLS Index 
# historical data series from Barclays website
# http://www.barclayhedge.com/cgi-bin/barclay_stats/bcndx.cgi?dump=excel&daily=daily&prog_cod=morningSta

# Peter Carl

# DETAILS
# Download and parse the Morningstar Long/Short Commodity Index

# Load needed packages:
require(zoo)
require(gdata)
require(quantmod)

# Download the first sheet in the xls workbook directly from the web site:
x = read.xls("http://www.barclayhedge.com/cgi-bin/barclay_stats/bcndx.cgi?dump=excel&daily=daily&prog_cod=morningSta", pattern="Date")
# x = read.xls("/home/peter/Data/Commodity Returns.xls", sheet="Daily", pattern="Date")

# That should result in something like this:
# > head(x)
#          Date     ROR    VAMI
# 1 Dec-26-1979 1.1753% 1011.75
# 2 Dec-27-1979 0.5799% 1017.62
# 3 Dec-28-1979 0.5103% 1022.81
# 4 Dec-31-1979 0.7056% 1030.03
# 5 Jan-01-1980 0.0000% 1030.03
# 6 Jan-02-1980 0.9884% 1040.21

# The first column is the date, in '%b-%d-%Y' format:
ISOdates = as.Date(x[,1], "%b-%d-%Y")
# ISOdates = as.Date(x[,"Date"], "%F") # new extract
# The second column is the return, already calculated:
ROR = as.numeric(gsub("[%]", "", x[,2]))/100

# Create a return series
x.R = xts(ROR, order.by=ISOdates)
colnames(x.R)= "Returns"
# Create a price series or index
x.IDX = xts(x[,3], order.by=ISOdates)
# x.IDX = xts(x[,"Index.Level"], order.by=ISOdates) # new extract
colnames(x.IDX) = "Close"
# new extract
x.R = Return.calculate(x.IDX)
colnames(x.R)= "Returns"
# Merge the series into a single object
MstarCLS.IDX = cbind(x.IDX, x.R)

# Calculate monthly returns from daily index values
x.M.IDX = to.monthly(Cl(MstarCLS.IDX))
colnames(x.M.IDX)= c("Open", "High", "Low", "Close")
x.M.R = Return.calculate(Cl(x.M.IDX))
colnames(x.M.R)="Returns"
MstarCLS.M.IDX = cbind(x.M.IDX, x.M.R)
index(MstarCLS.M.IDX) <- as.Date(index(MstarCLS.M.IDX), frac=1)

# Draw some pictures
charts.PerformanceSummary(MstarCLS.IDX["2010::","Returns"], ylog=TRUE, wealth.index=TRUE, main = "Morningstar CLS Index Returns")
MstarCLS.table=table.CalendarReturns(MstarCLS.M.IDX[,"Returns"])

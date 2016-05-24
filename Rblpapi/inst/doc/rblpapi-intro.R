## ---- eval=FALSE---------------------------------------------------------
#  library(Rblpapi)

## ---- eval=FALSE---------------------------------------------------------
#  blpConnect()

## ---- eval=FALSE---------------------------------------------------------
#  bdp(c("ESA Index", "SPY US Equity"), c("PX_LAST", "VOLUME"))

## ---- eval=FALSE---------------------------------------------------------
#  bds("GOOG US Equity", "TOP_20_HOLDERS_PUBLIC_FILINGS")

## ---- eval=FALSE---------------------------------------------------------
#  bdh("SPY US Equity", c("PX_LAST", "VOLUME"), start.date=Sys.Date()-31)

## ---- eval=FALSE---------------------------------------------------------
#  getBars("ES1 Index")

## ---- eval=FALSE---------------------------------------------------------
#  getTicks("ES1 Index")

## ---- eval=FALSE---------------------------------------------------------
#  res <- fieldSearch("vwap")

## ---- eval=FALSE---------------------------------------------------------
#  beqs("Global Oil Companies YTD Return","GLOBAL")


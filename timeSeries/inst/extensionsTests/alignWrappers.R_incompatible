

# artificial 1 sec data with missing Data


tX <- timeSequence("2014-03-07 00:00:00", "2014-03-07 23:59:59", by="sec")

s <- sample(1:length(tX))[1:length(tX)/10]
tX <- tX[-s]



###############################################################################
# align
#   extract index values of a given xts object corresponding to the last 
#   observations given a period specified by on


require(timeSeries)


# Random Seed:
set.seed(1953)

# Create a day of 1s time stamps:
tX <- timeSequence("2014-03-07 09:03:17", "2014-03-07 15:53:16", by="sec")

# Remove randomly 10% of the data:
s <- sample(1:length(tX))[1:length(tX)/10]
tX <- sort(tX[-s])
tS <- 201.7*cumulated(timeSeries(data=rnorm(length(tX))/(24*3600), charvec=tX))

plot(tS)
head(tS)


tZ <- align(tS, by="1min", method="fillNA", offset="42s")  
head(tZ)

tZ <- align(tS, by="3min", method="fillNA", offset="162s")  
head(tZ)

tZ <- align(tS, by="5min", method="fillNA", offset="102") 
head(tZ)

tZ <- align(tS, by="15min", method="fillNA", offset="702s")  
head(tZ)

tZ <- align(tS, by="30min", method="fillNA", offset="1602s")  
head(tZ)

tZ <- align(tS, by="60min", method="fillNA", offset="3402")  
head(tZ)


toPeriod <- function(x, by, method, offset="0s"")
{
  open <- function(x) as.vector(x)[1]
  high <- function(x) max(x)
  low <- function(x) min(x)
  close <- function(x) rev(as.vector(x))[1]

  cbind(
    aggregate(SPI, by, open),
    aggregate(SPI, by, high),
    aggregate(SPI, by, low),
    aggregate(SPI, by, close))
}

A1 <- timeSeries::align(tS, by="60min")
A2 <- xts::to.period(as.xts(tS), period = "minutes", k = 2) 


open <- function(x) as.vector(x)[1]
close <- function(x) rev(as.vector(x))[1]
high <- function(x) max(x)
low <- function(x) min(x)

SPI <- tS[, "SPI"]
by <- timeLastDayInMonth(time(tS))
OHLC <- cbind(
  aggregate(SPI, by, open),
  aggregate(SPI, by, high),
  aggregate(SPI, by, low),
  aggregate(SPI, by, close))
OHLC


xts::to.minutes(x,k,name,...)
xts::to.minutes3(x,name,...)
xts::to.minutes5(x,name,...)
xts::to.minutes10(x,name,...)
xts::to.minutes15(x,name,...)
xts::to.minutes30(x,name,...)
xts::to.hourly(x,name,...)


# -----------------------------------------------------------------------------

# Time alignment:

alignDaily(x=time(tS), include.weekends=FALSE)
alignMonthly(x=time(tS), include.weekends=FALSE)     # error
alignQuarterly(x=time(tS), include.weekends=FALSE)   # error


tD <- Sys.timeDate() + 1:1000
timeDate::align(tD, by="10s")
timeDate::align(tD, by="60s")
timeDate::align(tD, by="10m")     # error


td <- as.xts(Sys.time()) + 1:1000
xts::align.time(td, n=10)         # every 10 seconds
xts::align.time(td, n=60)         # align to next whole minute
xts::align.time(td, n=10*60)      # align to next whole 10 min interval


xts::shift.time(td, n=10)
xts::shift.time(td, n=60)

xts::shift.time(td)

# -----------------------------------------------------------------------------

xts::to.daily(x,drop.time=TRUE,name,...)

xts::to.weekly(x,drop.time=TRUE,name,...)
xts::to.monthly(x,indexAt='yearmon',drop.time=TRUE,name,...)
xts::to.quarterly(x,indexAt='yearqtr',drop.time=TRUE,name,...)
xts::to.yearly(x,drop.time=TRUE,name,...)

xts::to.period(
  x,
  period = 'months', 
  k = 1, 
  indexAt, 
  name=NULL,
  OHLC = TRUE,
  ...)


# -----------------------------------------------------------------------------


Convert an object to a specified periodicity lower than the given data 
object. For example, convert a daily series to a monthly series, or a 
monthly series to a yearly one, or a one minute series to an hourly 
series.


data(sample_matrix)
xts <- as.xts(sample_matrix)    # is daily

to.weekly(xts)
to.monthly(xts)
to.quarterly(xts)
to.yearly(xts)

tS <- as.timeSeries(sample_matrix)


% -----------------------------------------------------------------------------


as.numeric(as.POSIXct(time(tS)))
getFinCenter(tS)


indexTZ(xts, )
tzone(xts, )
tzone(xts) <- "GMT"
.index(xts, )


indexClass(xts)
class(time(tS))


% -----------------------------------------------------------------------------


.index <- function(x) as.numeric(as.POSIXct(time(x)))
.indexDate <- function(x) .index(x)%/%86400L
.indexday <- function(x) .index(x)%/%86400L
.indexmday <- function(x) as.POSIXlt(.POSIXct(.index(x)))$mday
.indexwday <- function(x) as.POSIXlt(.POSIXct(.index(x)))$wday
.indexweek <- function(x)
.indexmon <- function(x)
.indexyday <- function(x)
.indexyear <- function(x)

.indexhour <- function(x)
.indexmin <- function(x)
.indexsec <- function(x)


atoms

  
  
  


# Roll over fixed periods of length k point by point ...
# Functions borrowed from zoo
  
timeSeries::rollMin(
  x, k, na.pad = FALSE, align = c("center", "left", "right"),  ...) 
timeSeries::rollMax(
  x, k, na.pad = FALSE, align = c("center", "left", "right"),  ...)
timeSeries::rollMean(
  x, k, na.pad = FALSE, align = c("center", "left", "right"),  ...)
timeSeries::rollMedian(
  x, k, na.pad = FALSE, align = c("center", "left", "right"),  ...)
timeSeries::rollStats(
  x, k, FUN = mean, na.pad = FALSE, align = c("center", "left", "right"), ...) 


# Roll over Calendarical periods:

rollDailySeries(x, period="7d", FUN, ...)
rollMonthlySeries(x, period="12m", by="1m", FUN, ...)
# e.g. rollQuarterlySeries(x, period="12m", by="3m", FUN)
# e.g. rollYearlySeries
  
rollMonthlyWindows(x, period="12m", by="1m")
  
apply
applySeries


# period.apply
#   Apply a specified function to data over a given interval, where the 
#   interval is taken to be the data from INDEX[k] to INDEX[k+1], for 
#   k=1:(length(INDEX)-1).



x1 <- xts(matrix(1:(9*6),nc=6),
         order.by=as.Date(13000,origin="1970-01-01")+1:9)
x2 <- x1

xtsAttributes(x1) <- list(series1="1")
xtsAttributes(x2) <- list(series2="2")

xtsAttributes(x1)
xtsAttributes(x2)


x3 <- x1+x2
xtsAttributes(x3)

x33 <- cbind(x1, x2)
xtsAttributes(x33)

x33 <- rbind(x2, x1)
xtsAttributes(x33)



###############################################################################

appendList <- function (x, value) {
  stopifnot(is.list(x), is.list(value))
  xnames <- names(x)
  for (v in names(value)) {
    x[[v]] <- 
      if (v %in% xnames && is.list(x[[v]]) && is.list(value[[v]])) 
        appendList(x[[v]], value[[v]])
      else c(x[[v]], value[[v]]) }
  x }


"setAttributes<-" <- function(obj, value) {
  stopifnot(is.list(value))
  ATTRIBUTES <- getAttributes(obj)
  VALUE <- appendList(ATTRIBUTES, value)
  attr(obj@documentation, "Attributes") <- VALUE
  obj }


getAttributes <- function(obj) {
  attr(obj@documentation, "Attributes") }




obj1 <- dummySeries()
getAttributes(obj1)
setAttributes(obj1) <- list(series="obj1")
getAttributes(obj1)


obj2 <- dummySeries()
getAttributes(obj2)
setAttributes(obj2) <- list(series="obj2")
getAttributes(obj2)


getAttributes(obj1+obj2)  # returns the attributes only for the first
getAttributes(obj1-obj2)  # returns the attributes only for the first

getAttributes(cbind(obj1, obj2))
getAttributes(cbind(obj1, as.matrix(obj2)))   # matrix fails


getAttributes(rbind(obj1, obj2))
getAttributes(rbind(obj1, as.matrix(obj2)))   # matrix fails


getAttributes( rev(obj) )

getAttributes( obj[, 1] )

getAttributes( sample(obj) )

getAttributes( sort(sample(obj)) ) 

getAttributes( scale(obj) ) 




getAttributes( returns(obj) ) 
getAttributes( cumulated(returns(obj)) ) 










BIND(# Add another Attribute:
ATTRIBUTES <- attr(obj@documentation, "Attributes")
ATTRIBUTES

ATTRIBUTES <- appendList(ATTRIBUTES, list(say="hello"))
ATTRIBUTES

attr(obj@documentation, "Attributes") <- ATTRIBUTES

cbind(obj, obj, documentation = obj@documentation)





# Documentation

# Series:
#    dim(@.Data)
#    @units
#    @positions
#    @format
#    @FinCenter
#    @recordIDs
#    @title
#    @documentation
#       attributes(@documentation, "attributes)
    
    










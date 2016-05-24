

require(timeSeries)

X <- cumulated(LPP2005REC)[, 1:3]
for (i in 1:3) X[, i] <- 100*X[, i]/as.vector(X[1,i])
Data <- alignDailySeries(X)  # add: startDate
Index <- time(Data)

# Generate time Series:
tS <- timeSeries(data=Data, charvec=format(Index))
tR <- returns(tS)


###############################################################################
# aggregate:
#   from stats Package: The function aggregate splits the data into subsets, 
#   computes summary statistics for each, and returns the result in a 
#   convenient form.

# AGGREGATION OVER NON-OVEWRLAPPING PERIODS
#    starting point: aligned daily Data


# Aggregation Function:
#   function (x, by, FUN, ...) 


# Aggregation Levels:
#   weekly/biweekly:    endOfWeek, onTuesdays, lastBusinessDay
#   monthly:            endOMonth, lastFriday, lastBusinessDay
#   quarterly:          3-monthly
#   half-annually:      6-monthly
#   yearly:             12-monthly


# Aggregation Statistics:
#   mean, sd, var, median, ...
open <- function(x) as.vector(x)[1]
close <- function(x) rev(as.vector(x))[1]
high <- function(x) max(x)
low <- function(x) min(x)
spread <- function(x) max(x) - min(x)


# -----------------------------------------------------------------------------
# Weekly


# End-of-week:
by <- timeLastDayInMonth(time(tS))
mean.tR <- aggregate(tR[, "SPI"], by, mean)
sd.tR <- aggregate(tR[, "SPI"], by, sd)
plot(cbind(mean.tR, sd.tR), type="h")


# Weekly - Last Zurich Business Day In Week
by <- timeLastBizdayInMonth(time(tS), holidays = holidayZURICH())
mean.tR <- aggregate(tR[, "SPI"], by, mean)
sd.tR <- aggregate(tR[, "SPI"], by, sd)
plot(cbind(mean.tR, sd.tR), type="h")


# Weekly on Tuesdays
by <- timeSequence(from=start(tD), to=end(tD), by = "week")
mean.tR <- aggregate(tR[, "SPI"], by, mean)
sd.tR <- aggregate(tR[, "SPI"], by, sd)
plot(cbind(mean.tR, sd.tR), type="h")


# -----------------------------------------------------------------------------
# Monthly


# -----------------------------------------------------------------------------
# Quarterly



###############################################################################



# -----------------------------------------------------------------------------
# End-of-Month Statistics


# -----------------------------------------------------------------------------
# Monthly Open-High-Low-Close

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


# -----------------------------------------------------------------------------
# Monthly Spread / Percentual Spread

spread <- function(x) max(x) - min(x)
pspread <- function(x) (max(x) - min(x)) / (0.5 * (max(x) + min(x)))

SPI <- tS[, "SPI"]
by <- timeLastDayInMonth(time(tS))
SPREAD <- cbind(
  Points=aggregate(SPI, by, spread),
  Percent=100*aggregate(SPI, by, pspread))
SPREAD <- round(SPREAD, 2)
SPREAD


################################################################################
# Rolling: Aggregation with Overlappinng Periods


# -----------------------------------------------------------------------------
# rolling 52-weekly-highs and lows





# xts: Mean on weekly Periods
ep <- xts::endpoints(x.xts, on='weeks', k=1)
by1 <- index(x.xts)[ep[-1]]
period1 <- xts::period.apply(x.xts, INDEX=ep, FUN=mean)


  
###############################################################################
# xts::apply.monthly
  
FUN <- mean
x <- x.xts 
  
apply.daily(x, FUN)
apply.weekly(x, FUN)
apply.monthly(x, FUN)
apply.quarterly(x, FUN)
apply.yearly(x, FUN)
  
  
# timeDate::align

FUN <- mean
x <- unique(time(x.tS))

alignDaily(x, include.weekends=FALSE)
  
by1 <- unique(alignMonthly(x, include.weekends=FALSE))
x1 <- timeSeries::aggregate(x.tS, by1, FUN)
  
by2 <- unique(alignMonthly(x, include.weekends=TRUE))
x2 <- timeSeries::aggregate(x.tS, by2, FUN)
  
  
  
by1 <- unique(alignQuarterly(x, include.weekends=FALSE))
x1 <- timeSeries::aggregate(x.tS, by1, FUN)
  
by2 <- unique(alignQuarterly(x, include.weekends=TRUE))
x2 <- timeSeries::aggregate(x.tS, by2, FUN)

cbind(x1,x2)





############################################################################### 



xts::first(x.xts)
xts::last(x.xts)

first2 <- function(x) x[start(x), ]
last2 <- function(x) x[end(x), ]


first2(x.tS)
last2(x.tS)


# -----------------------------------------------------------------------------


INDEX <- seq(1, nrow(xts), by=21)
INDEX

.period.apply(tS, INDEX, FUN=max)

.period.max <- function(x, INDEX, FUN=max) .period.apply(x, INDEX, max)
.period.max(tS[, 1], INDEX)

.period.min <- function(x, INDEX) .period.apply(x, INDEX, min)
.period.min(tS[, 1], INDEX)


xts::period.apply(xts[, 1], INDEX, FUN=max)
xts::period.max(xts[, 1], INDEX)
xts::period.min(xts[, 1], INDEX)
xts::period.prod(xts[, 1], INDEX)
xts::period.sum(xts[, 1], INDEX)


# -----------------------------------------------------------------------------
# timeBased


is.timeBased <- 
  function (x) 
{
  if (!any(sapply(c(
    "Date", "POSIXt", "chron", "dates", "times", 
    "timeDate", "yearmon", "yearqtr", "xtime"), 
    function(xx) inherits(x, xx)))) 
  {
    ans <- FALSE
  } else {
    ans <- TRUE
  }
  ans
}


timeBased <- function(x) { is.timeBased(x) }


# -----------------------------------------------------------------------------


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


xts::to.minutes(x,k,name,...)
xts::to.minutes3(x,name,...)
xts::to.minutes5(x,name,...)
xts::to.minutes10(x,name,...)
xts::to.minutes15(x,name,...)
xts::to.minutes30(x,name,...)
xts::to.hourly(x,name,...)
xts::to.daily(x,drop.time=TRUE,name,...)

xts::to.weekly(x, drop.time=TRUE, name,...)
xts::to.monthly(x, indexAt='yearmon', drop.time=TRUE,name,...)
xts::to.quarterly(x, indexAt='yearqtr', drop.time=TRUE,name,...)
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


#Convert an object to a specified periodicity lower than the given data 
#object. For example, convert a daily series to a monthly series, or a 
#monthly series to a yearly one, or a one minute series to an hourly 
#series.


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



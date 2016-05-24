# Chapter 6 - Going on a Date with R

# Working with Dates

xd <- as.Date("2012-07-27")
xd
str(xd)
weekdays(xd)
xd + 7
xd + 0:6
weekdays(xd + 0:6)

startDate <- as.Date("2012-01-01")
xm <- seq(startDate, by="2 months", length.out=6)
xm
 months(xm)
quarters(xm)

Sys.localeconv()

as.Date("27 July 2012", format="%d %B %Y")

as.Date("27/7/12", format="%d/%m/%y")

# Adding Time Information to Dates

apollo <- "July 20, 1969, 20:17:39"
apollo.fmt <- "%B %d, %Y, %H:%M:%S"
xct <- as.POSIXct(apollo, format=apollo.fmt, tz="UTC")
xct

format(xct, "%d/%m/%y")
format(xct, "%S minutes past %I %p, on %d %B %Y")

# Performing Operations on Dates and Times

## Addition and subtraction

24*60*60
xct + 7*86400
xct + 3*60*60
xct - 7*86400
as.Date(xct) - 7

## Comparison of dates

Sys.time()
Sys.time() < xct

dec.start <- as.POSIXct("1950-01-01")
dec <- seq(dec.start, by="10 years", length.out=4)
dec
dec > xct

## Extraction

xlt <- as.POSIXlt(xct)
xlt
xlt$year
xlt$mon
unclass(xlt)










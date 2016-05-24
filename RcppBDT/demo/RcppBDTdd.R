
## demo for date duration objects

require(utils, quiet=TRUE, warn=FALSE)
require(Rcpp, quiet=TRUE, warn=FALSE)
require(RcppBDT, quiet=TRUE, warn=FALSE)

cat("\n(1) Create a bdtDd object and print it -- as a difftime()\n")
dd <- new(bdtDd, 5)
print(dd)  ## converts to difftime

cat("\n(2) Simple arithmetic, adding btdDd object and adding/subtracing integers\n")
print(dd + dd)                          # date_duration + date_duration
print(dd + 3)
print(dd - 2)
print(2L + dd)
print(2 + dd)
try(print(2 - dd))				# not permitted, bails inside the method -- should we test for ops?

cat("\n(3) Creating a bdtDt date object and adding the bdtDd object\n")
dt <- new(bdtDt, 2012, 10, 10)
print(dt)
print(dd + dt)
print(dt + dd)
print(dt - dd)

cat("\n(4) Adding to the bdtDt date object\n")
print(dt + dd + days(3) + weeks(2))
print(dt)
dt$addMonths(2)
dt$addYears(1)
print(dt)

## no converter from Date yet:  today <- Sys.Date()    # base R, provides Date
##print(today + dd)

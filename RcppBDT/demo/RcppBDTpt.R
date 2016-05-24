
## demo for time objects

require(utils, quiet=TRUE, warn=FALSE)
require(Rcpp, quiet=TRUE, warn=FALSE)
require(RcppBDT, quiet=TRUE, warn=FALSE)

op <- options("digits.secs"=6)          # to make sure R display microseconds

## use one of several bdtPt ctors
pt <- new(bdtPt, ISOdatetime(2012,1,1,00,30,0,tz="UTC"))

du <- new(bdtDu, 0, 1, 2, 123456789)

print(pt + du)
print(du + pt)

print(pt + hours(3) + minutes(4))
print(pt + 5.5)
print(pt + 5.005)
print(pt + 5.005005)
print(pt + 5.005005005)
print(5.5 + pt)

## Need to make this work too -- second argument is ignored right now
print(pt + c(5.5, 5.5005))

options(op)


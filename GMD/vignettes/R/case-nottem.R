## load library
require("GMD")

## load data
data(nottem)

class(nottem)      # a time-series (ts) object
x <- ts2df(nottem) # convert ts to data.frame
mhist1 <- as.mhist(x[1:3,])

## plot multiple discrete distributions side-by-side
plot(mhist1,xlab="Month",ylab="Degree Fahrenheit",
     main="Air temperatures at Nottingham Castle")

## make summary statistics for each bin
mhist2 <- as.mhist(x)
ms <- mhist.summary(mhist2)
print(ms)

## plot bin-wise summary statistics with
## confidence intervals over the bars
plot(ms, main="Mean air temperatures at Nottingham Castle (1920-1939)",
     xlab="Month", ylab="Degree Fahrenheit")

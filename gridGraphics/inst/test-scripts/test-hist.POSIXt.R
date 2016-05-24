
library(gridGraphics)

hist.POSIXt1 <- function() {
    hist(.leap.seconds, "years", freq = TRUE)
}

hist.POSIXt2 <- function() {
    hist(.leap.seconds[.leap.seconds < ISOdate(2020, 1, 1)],
         seq(ISOdate(1970, 1, 1), ISOdate(2020, 1, 1), "5 years"))
}

hist.POSIXt3 <- function() {
    ## 100 random dates in a 10-week period
    set.seed(1)
    random.dates <- as.Date("2001/1/1") + 70*stats::runif(100)
    hist(random.dates, "weeks", format = "%d %b")
}

plotdiff(expression(hist.POSIXt1()), "hist.POSIXt-1", width=15)
plotdiff(expression(hist.POSIXt2()), "hist.POSIXt-2")
plotdiff(expression(hist.POSIXt3()), "hist.POSIXt-3", width=10)

plotdiffResult()

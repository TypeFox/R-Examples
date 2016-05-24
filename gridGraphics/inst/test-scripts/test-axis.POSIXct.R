
library(gridGraphics)

attach(beaver1)

axis.POSIXct1 <- function() {
    time <- strptime(paste(1990, day, time %/% 100, time %% 100),
                     "%Y %j %H %M")
    plot(time, temp, type = "l") # axis at 4-hour intervals.
}

axis.POSIXct2 <- function() {
    time <- strptime(paste(1990, day, time %/% 100, time %% 100),
                     "%Y %j %H %M")
    # now label every hour on the time axis
    plot(time, temp, type = "l", xaxt = "n")
    r <- as.POSIXct(round(range(time), "hours"))
    axis.POSIXct(1, at = seq(r[1], r[2], by = "hour"), format = "%H")
}

axis.POSIXct3 <- function() {
    plot(.leap.seconds, seq_along(.leap.seconds), type = "n", yaxt = "n",
         xlab = "leap seconds", ylab = "", bty = "n")
    rug(.leap.seconds)
}

axis.POSIXct4 <- function() {
    ## or as dates
    lps <- as.Date(.leap.seconds)
    plot(lps, seq_along(.leap.seconds),
         type = "n", yaxt = "n", xlab = "leap seconds",
         ylab = "", bty = "n")
    rug(lps)
}

axis.POSIXct5 <- function() {
    ## 100 random dates in a 10-week period
    set.seed(1)
    random.dates <- as.Date("2001/1/1") + 70*sort(stats::runif(100))
    plot(random.dates, 1:100)
}

axis.POSIXct6 <- function() {
    # or for a better axis labelling
    set.seed(1)
    random.dates <- as.Date("2001/1/1") + 70*sort(stats::runif(100))
    plot(random.dates, 1:100, xaxt = "n")
    axis.Date(1, at = seq(as.Date("2001/1/1"), max(random.dates)+6, "weeks"))
    axis.Date(1, at = seq(as.Date("2001/1/1"), max(random.dates)+6, "days"),
              labels = FALSE, tcl = -0.2)
}

plotdiff(expression(axis.POSIXct1()), "axis.POSIXct-1")
plotdiff(expression(axis.POSIXct2()), "axis.POSIXct-2", width=8, height=8)
plotdiff(expression(axis.POSIXct3()), "axis.POSIXct-3")
plotdiff(expression(axis.POSIXct4()), "axis.POSIXct-4")
plotdiff(expression(axis.POSIXct5()), "axis.POSIXct-5")
plotdiff(expression(axis.POSIXct6()), "axis.POSIXct-6", width=9, height=9)

detach(beaver1)

plotdiffResult()


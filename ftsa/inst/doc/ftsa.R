### R code from vignette source 'ftsa.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: ftsa.Rnw:112-115
###################################################
library(ftsa)
plot(forecast(ftsm(Australiasmoothfertility,
    order=2), h=20), "components")


###################################################
### code chunk number 2: ftsa.Rnw:144-151
###################################################
# Plot the historical data in gray
plot(Australiasmoothfertility, col = gray(0.8), xlab = "Age",
    ylab = "Count of live birth (per 1,000 females)",
    main = "Forecasted fertility rates (2007-2026)")
# Plot the forecasts in rainbow color for Figure 4(a)
plot(forecast(ftsm(Australiasmoothfertility, order=2), h = 20), add = TRUE)
legend("topright", c("2007", "2026"), col = c("red", "blue"), lty = 1)


###################################################
### code chunk number 3: ftsa.Rnw:181-189
###################################################
# Plot the point forecast
aus = forecast(ftsm(Australiasmoothfertility, order=2), h=1)
plot(aus, ylim=c(0,140))
# Plot the lower and upper bounds
lines(aus$lower, col=2)
lines(aus$upper, col=2)
# Add a legend to the plot
legend("topright", c("Point forecasts", "Interval forecasts"),col=c(1,2), lty=1, cex=0.9)


###################################################
### code chunk number 4: ftsa.Rnw:263-277
###################################################
# Name history to represent historical data,
history = ElNino2011smooth
# Name obs to represent partially observed data,
obs = ElNino2011smooth$y[1:5,62]
# Name fore to represent the forecasting period
fore = ElNino2011smooth$y[6:12,62]
int = dynupdate(data=history, newdata=obs, holdoutdata=fore, method="block", interval=TRUE)
bmupdate = dynupdate(data=history, newdata=obs, holdoutdata=fore, method="block", value=TRUE)
plot(6:12, fore, type="l", ylim=c(19,26),  xlab="Month", ylab="Sea surface temperature")
lines(6:12, bmupdate, col=4)
lines(6:12, int$low$y, col=2)
lines(6:12, int$up$y, col=2)
legend("topright", c("True observations", "Point forecasts", "Interval forecasts"), 
       col=c(1,4,2), lty=1, cex=0.8)



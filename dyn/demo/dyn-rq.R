
library(dyn)
library(quantreg)

# create zoo series from airquality Ozone column
ozone.z <- zoo(airquality[,"Ozone"])

# quantile regression of ozone on one day lag
ozone.rq <- dyn$rq(ozone.z ~ lag(ozone.z,-1), 0:10/10)

# plot data and overlay quantile regression lines
plot(formula(ozone.rq), pch = 20)
apply(coef(ozone.rq), 2, abline)

# plot each coefficient vs. tau, one graph per coefficient
plot(zoo(t(coef(ozone.rq)), ozone.rq$tau), xlab = "tau", main = "ozone")



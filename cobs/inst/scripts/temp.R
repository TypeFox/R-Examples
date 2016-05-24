#### Example 1 of   He and Ng (1999) :
#### -------------------------------
if(!interactive()) postscript("temp.ps", horizontal = TRUE)

library(cobs)

data(globtemp)
str(globtemp)## time series T=113,  1880:1992

## require(ts)# later for plotting
## At the moment forget about all time-series stuff; use num.vectors:
str(year <- as.vector(time(globtemp)))
str(temp <- as.vector(globtemp))

## Codes for Figure 1a
a50 <- cobs(year, temp, knots.add = TRUE, degree = 1, constraint = "increase")
## As suggested in the warning message, we increase the number of knots to 9
a50 <- cobs(year, temp, nknots = 9, knots.add = TRUE, degree = 1,
            constraint = "increase")

## Here, we use the same knots sequence chosen for the 50th percentile
a10 <- cobs(year, temp, nknots = length(a50$knots), knots = a50$knot,
            degree = 1, tau = 0.1, constraint = "increase")
a90 <- cobs(year, temp, nknots = length(a50$knots), knots = a50$knot,
            degree = 1, tau = 0.9, constraint = "increase")
hot.idx  <- temp >= a90$fit
cold.idx <- temp <= a10$fit
normal.idx <- !hot.idx & !cold.idx

pa50 <- predict(a50,year,interval = "none")
pa10 <- predict(a10,year,interval = "none")
pa90 <- predict(a90,year,interval = "none")

plot(year, temp, type = "n", ylab = "Temperature (C)", ylim = c(-.7,.6))
lines(pa50, col = 2)
lines(pa10, col = 3)
lines(pa90, col = 3)
points(year, temp, pch = c(1,3)[2 - normal.idx])

text(year[hot.idx], temp[hot.idx] + .03, labels = year[hot.idx])
text(year[cold.idx],temp[cold.idx]- .03, labels = year[cold.idx])

### Codes for Figure 1b. -- UNconstrained

## (Pin Ng. had these commented out)

a50 <- cobs(year, temp, knots.add = TRUE, degree = 1, constraint = "none")
## As suggested in the warning message, we increase the number of knots to 9
a50 <- cobs(year, temp, nknots = 9, knots.add = TRUE, degree = 1,
            constraint = "none")
## Here, we use the same knots sequence chosen for the 50th percentile
a10 <- cobs(year, temp, nknots = length(a50$knots), knots = a50$knot,
            degree = 1, tau = 0.1, constraint = "none")
a90 <- cobs(year, temp, nknots = length(a50$knots), knots = a50$knot,
            degree = 1, tau = 0.9, constraint = "none")
hot.idx <- temp >= a90$fit
cold.idx <- temp <= a10$fit
normal.idx <- (temp < a90$fit ) & (temp > a10$fit)
pa50 <- predict(a50,year,interval = "none")
pa10 <- predict(a10,year,interval = "none")
pa90 <- predict(a90,year,interval = "none")
plot(year, temp, type = "n",xlab = "Year", ylab = "Temperature (C)",
     ylim = c(-.7,.6))
lines(pa50[,1], pa50[,2])
lines(pa10[,1], pa10[,2])
lines(pa90[,1], pa90[,2])
points(year[normal.idx],temp[normal.idx])
points(year[!normal.idx],temp[!normal.idx],pch = 3)
text(year[hot.idx],temp[hot.idx]+.05,labels = year[hot.idx])
text(year[cold.idx],temp[cold.idx]-.05,labels = year[cold.idx])


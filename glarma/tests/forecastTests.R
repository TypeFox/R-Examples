require(glarma)
require(zoo)
### Model number of deaths
data(DriverDeaths)
y <- DriverDeaths[, "Deaths"]
X <- as.matrix(DriverDeaths[, 2:5])
Population <- DriverDeaths[, "Population"]

### Offset included
glarmamod <- glarma(y, X, offset = log(Population/100000),
                    phiLags = c(12),
                    thetaLags = c(1),
                    type = "Poi", method = "FS",
                    residuals = "Pearson", maxit = 100, grad = 1e-6)
print(summary(glarmamod))

XT1 <- matrix(X[72,], nrow = 1)
offsetT1 <- log(Population/100000)[72]

mu <- forecast(glarmamod, 1, XT1, offsetT1)$mu
print(mu)


### Save some values
allX <- X
allFits <- fitted(glarmamod)
ally <- y

### Look at a succession of forecasts
### Using actual values in forecasts
forecasts <- numeric(72)
for (i in (62:71)){
    y <- DriverDeaths[1:i, "Deaths"]
    X <- as.matrix(DriverDeaths[1:i, 2:5])
    Population <- DriverDeaths[1:i, "Population"]

    ## Offset included
    glarmamod <- glarma(y, X, offset = log(Population/100000),
                        phiLags = c(12),
                        thetaLags = c(1),
                        type = "Poi", method = "FS",
                        residuals = "Pearson", maxit = 100, grad = 1e-6)
    XT1 <- matrix(allX[i + 1, ], nrow = 1)
    offsetT1 <- log(DriverDeaths$Population[i + 1]/100000)
    mu <- forecast(glarmamod, 1, XT1, offsetT1)$mu
    if (i == 62){
        forecasts[1:62] <- fitted(glarmamod)
    }
    forecasts[i+1] <- mu
}
par(mfrow = c(1,1))
forecasts <- ts(forecasts[63:72], start = c(1985, 10), deltat = 1/12)
fitted <- ts(allFits, start = c(1980, 8), deltat = 1/12)
obs <- ts(DriverDeaths$Deaths, start = c(1980, 8), deltat = 1/12)
plot(obs, ylab = "Driver Deaths", lty = 2,
     main = "Single Vehicle Nighttime Driver Deaths in Utah")
points(obs)
lines(fitted, lwd = 2)
lines(forecasts, col = "red")
par(xpd = NA)
graph.param <-
    legend("top",
           legend = c("observations",expression(estimated~mu[t]),
                      expression(predicted~mu[t])),
           ncol = 3,
           cex = 0.7,
           bty = "n", plot = FALSE)
legend(graph.param$rect$left,
       graph.param$rect$top + graph.param$rect$h,
       legend = c("observations", expression(estimated~mu[t]),
                  expression(predicted~mu[t])),
       col = c("black","black","red"),
       lwd = c(1,2,1), lty = c(2,1,1),
       pch = c(1, NA_integer_, NA_integer_),
       ncol = 3,
       cex = 0.7,
       bty = "n",
       text.font = 4)
par(xpd = FALSE)

### Polio:
data(Polio)
y  <-  Polio[, 2]
X  <-  as.matrix(Polio[, 3:8])
glarmamod <- glarma(y, X, thetaLags = c(1,2,5),
                    type = "NegBin", method = "FS",
                    residuals = "Score", maxit = 100, grad = 1e-6)

summary(glarmamod)

XT1 <- matrix(X[168, ], nrow = 1)
mu <- forecast(glarmamod, 1, XT1, offsetT1 = 0)$mu


plot(0:11, dnbinom(0:11, mu = mu, size =  coef(glarmamod)$NB))

### Save some values
allX <- X
allFits <- fitted(glarmamod)
ally <- y

### Look at a succession of forecasts
### Using actual values in forecasts
forecasts <- numeric(168)
for (i in (156:167)){
    y  <-  Polio[1:i, 2]
    X  <-  as.matrix(Polio[1:i, 3:8])
    glarmamod <- glarma(y, X, thetaLags = c(1,2,5),
                        type = "NegBin", method = "FS",
                        residuals = "Score", maxit = 100, grad = 1e-6)
    XT1 <- matrix(c(1, (i + 1 - 72)/1000, cos(2*pi*(i + 1)/12),
             sin(2*pi*(i + 1)/12), cos(4*pi*(i + 1)/12),
             sin(4*pi*(i + 1)/12)), nrow = 1)
    offsetT1 <- 0
    mu <- forecast(glarmamod, 1, XT1, offsetT1)$mu
    if (i == 156){
        forecasts[1:156] <- fitted(glarmamod)
    }
    forecasts[i+1] <- mu
}
par(mfrow = c(1,1))
forecasts <- ts(forecasts[157:168], start = c(1983, 1), deltat = 1/12)
fitted <- ts(allFits, start = c(1970, 1), deltat = 1/12)
data <- ts(Polio$Cases, start = c(1970, 1), deltat = 1/12)
plot(data, ylab = "Polio Cases",
     main = "Monthly Numbers of Polio Cases in the USA")
points(data)
lines(fitted, lwd = 2)
lines(forecasts, col = "red")
par(xpd = NA)
graph.param <-
    legend("top",
           legend = c("observations",expression(estimated~mu[t]),
                      expression(predicted~mu[t])),
           lty = c(2,1,1),
           ncol = 3,
           cex = 0.7,
           bty = "n", plot = FALSE)
legend(graph.param$rect$left,
       graph.param$rect$top + graph.param$rect$h,
       legend = c("observations", expression(estimated~mu[t]),
                  expression(predicted~mu[t])),
       col = c("black","black","red"),
       lwd = c(1,2,1), lty = c(2,1,1),
       pch = c(1, NA_integer_, NA_integer_),
       ncol = 3,
       cex = 0.7,
       bty = "n",
       text.font = 4)
par(xpd = FALSE)

### Boat Race:
data(OxBoatRace)

y1 <- OxBoatRace$Camwin
n1 <- rep(1, length(OxBoatRace$Year))
Y <- cbind(y1, n1 - y1)
X <- cbind(OxBoatRace$Intercept, OxBoatRace$Diff)
colnames(X) <- c("Intercept", "Weight Diff")

oxcamglm <- glm(Y ~ Diff + I(Diff^2),
                data = OxBoatRace,
                family = binomial(link = "logit"), x = TRUE)

X <- oxcamglm$x

glarmamod <- glarma(Y, X, thetaLags = c(1, 2), type = "Bin", method = "NR",
                    residuals = "Pearson", maxit = 100, grad = 1e-6)

summary(glarmamod)

par(mfrow=c(3,2))
plot(glarmamod)


### Save some values
allX <- X
allFits <- fitted(glarmamod)
ally1 <- y1


### Look at a succession of forecasts
### Using actual values in forecasts
forecasts <- numeric(156)
for (i in (146:155)){
    y1 <- OxBoatRace$Camwin[1:i]
    n1 <- rep(1, i)
    Y <- cbind(y1, n1 - y1)
    X <- cbind(OxBoatRace$Intercept[1:i], OxBoatRace$Diff[1:i])
    colnames(X) <- c("Intercept", "Weight Diff")

    oxcamglm <- glm(Y ~ Diff + I(Diff^2),
                    data = OxBoatRace[1:i, ],
                    family = binomial(link = "logit"), x = TRUE)
    X <- oxcamglm$x
    glarmamod <- glarma(Y, X, thetaLags = c(1, 2), type = "Bin", method = "NR",
                    residuals = "Pearson", maxit = 100, grad = 1e-6)
    XT1 <- matrix(allX[i + 1, ], nrow = 1)
    offsetT1 <- 0
    mu <- forecast(glarmamod, 1, XT1, offsetT1 = 0)$mu
    if (i == 146){
        forecasts[1:146] <- fitted(glarmamod)
    }
    forecasts[i+1] <- mu
}

par(mfrow = c(1,1))
forecasts <- zoo(forecasts, order.by = OxBoatRace$Year)
fitted <- zoo(allFits, order.by = OxBoatRace$Year)
plot(fitted, lwd = 2, xlab = "Year", ylab = "Cambridge Win")
lines(OxBoatRace$Year[147:156], forecasts[147:156], col = "red")
plot(fitted[137:156], lwd = 2, ylim = c(0,1),
     xlab = "Year", ylab = "Cambridge Win",
     main = "Oxford versus Cambridge Boat Race")
lines(OxBoatRace$Year[147:156], forecasts[147:156], col = "red")
points(OxBoatRace$Year[137:156],OxBoatRace$Camwin[137:156])
abline(h=0.5)
par(xpd = NA)
graph.param <-
    legend("top",
           legend = c("observations",expression(estimated~mu[t]),
                      expression(predicted~mu[t])),
           lty = c(2,1,1),
           ncol = 3,
           cex = 0.7,
           bty = "n", plot = FALSE)
legend(graph.param$rect$left,
       graph.param$rect$top + graph.param$rect$h,
       legend = c("observations", expression(estimated~mu[t]),
                  expression(predicted~mu[t])),
       col = c("black","black","red"),
       lwd = c(NA,2,1), lty = c(NA,1,1),
       pch = c(1, NA_integer_, NA_integer_),
       ncol = 3,
       cex = 0.7,
       bty = "n",
       text.font = 4)
par(xpd = FALSE)



### is there skill in the GLARMA model?
results <- cbind(round(forecasts[147:156]),OxBoatRace$Camwin[147:156])
correct <- (results[,1] - results[,2] == 0)
mean(correct)

###############################
### Test multiple steps ahead
###############################
### Model number of deaths
data(DriverDeaths)
y <- DriverDeaths[, "Deaths"]
X <- as.matrix(DriverDeaths[, 2:5])
Population <- DriverDeaths[, "Population"]

### Offset included
glarmamod <- glarma(y, X, offset = log(Population/100000),
                    phiLags = c(12),
                    thetaLags = c(1),
                    type = "Poi", method = "FS",
                    residuals = "Pearson", maxit = 100, grad = 1e-6)
XT1 <- as.matrix(X[61:72, ])
offsetT1 <- log(Population/100000)[61:72]

nObs <- NROW(X)
n.ahead <- 10
XT1 <- as.matrix(X[(nObs - n.ahead +1):nObs, ])
offsetT1 <- log(Population/100000)[(nObs - n.ahead +1):nObs]
nSims <- 500
forecastY <- matrix(ncol = n.ahead, nrow = nSims)
forecastMu <- matrix(ncol = n.ahead, nrow = nSims)

for(i in 1:nSims){
    temp <-  forecast(glarmamod, n.ahead, XT1, offsetT1)
    forecastY[i, ] <- temp$Y
    forecastMu[i, ] <- temp$mu
}
par(mfrow = c(4, 3))
for (lead in 1:10) {
    hist(forecastY[, lead], freq = FALSE, breaks = 0:15,
         main = as.character(lead))
    print(table(forecastY[, lead]))
}

##########################
### Multiple steps ahead example for man page
### Generate a sample of Y values 2 steps ahead and examine the distribution
data(DriverDeaths)
y <- DriverDeaths[, "Deaths"]
X <- as.matrix(DriverDeaths[, 2:5])
Population <- DriverDeaths[, "Population"]

### Fit the glarma model to the first 70 observations
glarmamod <- glarma(y[1:70], X[1:70, ],
                    offset = log(Population/100000)[1:70],
                    phiLags = c(12),
                    thetaLags = c(1),
                    type = "Poi", method = "FS",
                    residuals = "Pearson", maxit = 100, grad = 1e-6)

nObs <- NROW(X)
n.ahead <- 2
### Specify the X matrix and offset for the times where predictions
### are required
XT1 <- as.matrix(X[(nObs - n.ahead + 1):nObs, ])
offsetT1 <- log(Population/100000)[(nObs - n.ahead + 1):nObs]
nSims <- 500
forecastY <- matrix(ncol = n.ahead, nrow = nSims)
forecastMu <- matrix(ncol = n.ahead, nrow = nSims)

### Generate sample predicted values
for(i in 1:nSims){
    temp <-  forecast(glarmamod, n.ahead, XT1, offsetT1)
    forecastY[i, ] <- temp$Y
    forecastMu[i, ] <- temp$mu
}
### Examine distribution of sample of Y values n.ahead
table(forecastY[, 2])
par(mfrow = c(2,1))
barplot(table(forecastY[, 2]),
        main = "Barplot of Sample Y Values 2 Steps Ahead")
hist(forecastY[, 2], xlab = "Sample Y values",
     main = "Histogram of Sample Y Values 2 Steps Ahead\nwith 0.025 and 0.975 Quantiles")
abline(v = quantile(forecastY[, 2], c(0.025, 0.975)), col = "red")

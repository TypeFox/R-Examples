###################################################
### chunk number 1: setup
###################################################
options(prompt = "R> ", continue = "+  ", width = 64,
  digits = 4, show.signif.stars = FALSE, useFancyQuotes = FALSE)

options(SweaveHooks = list(onefig =   function() {par(mfrow = c(1,1))},
                           twofig =   function() {par(mfrow = c(1,2))},                           
                           threefig = function() {par(mfrow = c(1,3))},
                           fourfig =  function() {par(mfrow = c(2,2))},
			   sixfig =   function() {par(mfrow = c(3,2))}))

library("AER")

set.seed(1071)


###################################################
### chunk number 2: options
###################################################
options(digits = 6)


###################################################
### chunk number 3: ts-plot eval=FALSE
###################################################
## data("UKNonDurables")
## plot(UKNonDurables)


###################################################
### chunk number 4: UKNonDurables-data
###################################################
data("UKNonDurables")


###################################################
### chunk number 5: tsp
###################################################
tsp(UKNonDurables)


###################################################
### chunk number 6: window
###################################################
window(UKNonDurables, end = c(1956, 4))


###################################################
### chunk number 7: filter eval=FALSE
###################################################
## data("UKDriverDeaths")
## plot(UKDriverDeaths)
## lines(filter(UKDriverDeaths, c(1/2, rep(1, 11), 1/2)/12),
##   col = 2)


###################################################
### chunk number 8: ts-plot1
###################################################
data("UKNonDurables")
plot(UKNonDurables)
data("UKDriverDeaths")
plot(UKDriverDeaths)
lines(filter(UKDriverDeaths, c(1/2, rep(1, 11), 1/2)/12),
  col = 2)


###################################################
### chunk number 9: filter1 eval=FALSE
###################################################
## data("UKDriverDeaths")
## plot(UKDriverDeaths)
## lines(filter(UKDriverDeaths, c(1/2, rep(1, 11), 1/2)/12),
##   col = 2)


###################################################
### chunk number 10: rollapply
###################################################
plot(rollapply(UKDriverDeaths, 12, sd))


###################################################
### chunk number 11: ar-sim
###################################################
set.seed(1234)
x <- filter(rnorm(100), 0.9, method = "recursive")


###################################################
### chunk number 12: decompose
###################################################
dd_dec <- decompose(log(UKDriverDeaths))
dd_stl <- stl(log(UKDriverDeaths), s.window = 13)


###################################################
### chunk number 13: decompose-components
###################################################
plot(dd_dec$trend, ylab = "trend")
lines(dd_stl$time.series[,"trend"], lty = 2, lwd = 2)


###################################################
### chunk number 14: seat-mean-sd
###################################################
plot(dd_dec$trend, ylab = "trend")
lines(dd_stl$time.series[,"trend"], lty = 2, lwd = 2)
plot(rollapply(UKDriverDeaths, 12, sd))


###################################################
### chunk number 15: stl
###################################################
plot(dd_stl)


###################################################
### chunk number 16: Holt-Winters
###################################################
dd_past <- window(UKDriverDeaths, end = c(1982, 12))
dd_hw <- HoltWinters(dd_past)
dd_pred <- predict(dd_hw, n.ahead = 24)


###################################################
### chunk number 17: Holt-Winters-plot
###################################################
plot(dd_hw, dd_pred, ylim = range(UKDriverDeaths))
lines(UKDriverDeaths)


###################################################
### chunk number 18: Holt-Winters-plot1
###################################################
plot(dd_hw, dd_pred, ylim = range(UKDriverDeaths))
lines(UKDriverDeaths)


###################################################
### chunk number 19: acf eval=FALSE
###################################################
## acf(x)
## pacf(x)


###################################################
### chunk number 20: acf1
###################################################
acf(x, ylim = c(-0.2, 1))
pacf(x, ylim = c(-0.2, 1))


###################################################
### chunk number 21: ar
###################################################
ar(x)


###################################################
### chunk number 22: window-non-durab
###################################################
nd <- window(log(UKNonDurables), end = c(1970, 4))


###################################################
### chunk number 23: non-durab-acf
###################################################
acf(diff(nd), ylim = c(-1, 1))
pacf(diff(nd), ylim = c(-1, 1))
acf(diff(diff(nd, 4)), ylim = c(-1, 1))
pacf(diff(diff(nd, 4)), ylim = c(-1, 1))


###################################################
### chunk number 24: non-durab-acf1
###################################################
acf(diff(nd), ylim = c(-1, 1))
pacf(diff(nd), ylim = c(-1, 1))
acf(diff(diff(nd, 4)), ylim = c(-1, 1))
pacf(diff(diff(nd, 4)), ylim = c(-1, 1))


###################################################
### chunk number 25: arima-setup
###################################################
nd_pars <- expand.grid(ar = 0:2, diff = 1, ma = 0:2,
  sar = 0:1, sdiff = 1, sma = 0:1)
nd_aic <- rep(0, nrow(nd_pars))
for(i in seq(along = nd_aic)) nd_aic[i] <- AIC(arima(nd,
  unlist(nd_pars[i, 1:3]), unlist(nd_pars[i, 4:6])),
  k = log(length(nd)))
nd_pars[which.min(nd_aic),]


###################################################
### chunk number 26: arima
###################################################
nd_arima <- arima(nd, order = c(0,1,1), seasonal = c(0,1,1))
nd_arima


###################################################
### chunk number 27: tsdiag
###################################################
tsdiag(nd_arima)


###################################################
### chunk number 28: tsdiag1
###################################################
tsdiag(nd_arima)


###################################################
### chunk number 29: arima-predict
###################################################
nd_pred <- predict(nd_arima, n.ahead = 18 * 4)


###################################################
### chunk number 30: arima-compare
###################################################
plot(log(UKNonDurables))
lines(nd_pred$pred, col = 2)


###################################################
### chunk number 31: arima-compare1
###################################################
plot(log(UKNonDurables))
lines(nd_pred$pred, col = 2)


###################################################
### chunk number 32: pepper
###################################################
data("PepperPrice")
plot(PepperPrice, plot.type = "single", col = 1:2)
legend("topleft", c("black", "white"), bty = "n", 
col = 1:2, lty = rep(1,2))


###################################################
### chunk number 33: pepper1
###################################################
data("PepperPrice")
plot(PepperPrice, plot.type = "single", col = 1:2)
legend("topleft", c("black", "white"), bty = "n", 
col = 1:2, lty = rep(1,2))


###################################################
### chunk number 34: adf1
###################################################
library("tseries")
adf.test(log(PepperPrice[, "white"]))


###################################################
### chunk number 35: adf1
###################################################
adf.test(diff(log(PepperPrice[, "white"])))


###################################################
### chunk number 36: pp
###################################################
pp.test(log(PepperPrice[, "white"]), type = "Z(t_alpha)")


###################################################
### chunk number 37: urca eval=FALSE
###################################################
## library("urca")
## pepper_ers <- ur.ers(log(PepperPrice[, "white"]),
##   type = "DF-GLS", model = "const", lag.max = 4)
## summary(pepper_ers)


###################################################
### chunk number 38: kpss
###################################################
kpss.test(log(PepperPrice[, "white"]))


###################################################
### chunk number 39: po
###################################################
po.test(log(PepperPrice))


###################################################
### chunk number 40: joh-trace
###################################################
library("urca")
pepper_jo <- ca.jo(log(PepperPrice), ecdet = "const",
  type = "trace")
summary(pepper_jo)


###################################################
### chunk number 41: joh-lmax eval=FALSE
###################################################
## pepper_jo2 <- ca.jo(log(PepperPrice), ecdet = "const", type = "eigen")
## summary(pepper_jo2)


###################################################
### chunk number 42: dynlm-by-hand
###################################################
dd <- log(UKDriverDeaths)
dd_dat <- ts.intersect(dd, dd1 = lag(dd, k = -1),
  dd12 = lag(dd, k = -12))
lm(dd ~ dd1 + dd12, data = dd_dat)


###################################################
### chunk number 43: dynlm
###################################################
library("dynlm")
dynlm(dd ~ L(dd) + L(dd, 12))


###################################################
### chunk number 44: efp
###################################################
library("strucchange")
dd_ocus <- efp(dd ~ dd1 + dd12, data = dd_dat,
  type = "OLS-CUSUM")


###################################################
### chunk number 45: efp-test
###################################################
sctest(dd_ocus)


###################################################
### chunk number 46: efp-plot eval=FALSE
###################################################
## plot(dd_ocus)


###################################################
### chunk number 47: Fstats
###################################################
dd_fs <- Fstats(dd ~ dd1 + dd12, data = dd_dat, from = 0.1)
plot(dd_fs)
sctest(dd_fs)


###################################################
### chunk number 48: ocus-supF
###################################################
plot(dd_ocus)
plot(dd_fs, main = "supF test")


###################################################
### chunk number 49: GermanM1
###################################################
data("GermanM1")
LTW <- dm ~ dy2 + dR + dR1 + dp + m1 + y1 + R1 + season


###################################################
### chunk number 50: re eval=FALSE
###################################################
## m1_re <- efp(LTW, data = GermanM1, type = "RE")
## plot(m1_re)


###################################################
### chunk number 51: re1
###################################################
m1_re <- efp(LTW, data = GermanM1, type = "RE")
plot(m1_re)


###################################################
### chunk number 52: dating
###################################################
dd_bp <- breakpoints(dd ~ dd1 + dd12, data = dd_dat, h = 0.1)


###################################################
### chunk number 53: dating-coef
###################################################
coef(dd_bp, breaks = 2)


###################################################
### chunk number 54: dating-plot eval=FALSE
###################################################
## plot(dd)
## lines(fitted(dd_bp, breaks = 2), col = 4)
## lines(confint(dd_bp, breaks = 2))


###################################################
### chunk number 55: dating-plot1
###################################################
plot(dd_bp, legend = FALSE, main = "")
plot(dd)
lines(fitted(dd_bp, breaks = 2), col = 4)
lines(confint(dd_bp, breaks = 2))


###################################################
### chunk number 56: StructTS
###################################################
dd_struct <- StructTS(log(UKDriverDeaths))


###################################################
### chunk number 57: StructTS-plot eval=FALSE
###################################################
## plot(cbind(fitted(dd_struct), residuals(dd_struct)))


###################################################
### chunk number 58: StructTS-plot1
###################################################
dd_struct_plot <- cbind(fitted(dd_struct), residuals = residuals(dd_struct))
colnames(dd_struct_plot) <- c("level", "slope", "season", "residuals")
plot(dd_struct_plot, main = "")


###################################################
### chunk number 59: garch-plot
###################################################
data("MarkPound")
plot(MarkPound, main = "")


###################################################
### chunk number 60: garch
###################################################
mp <- garch(MarkPound, grad = "numerical", trace = FALSE)
summary(mp)



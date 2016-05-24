### R code from vignette source 'CADFtest.Rnw'

###################################################
### code chunk number 1: CADFtest.Rnw:8-9
###################################################
options(prompt = "R> ", continue="+ ", useFancyQuotes=FALSE)


###################################################
### code chunk number 2: CADFtest.Rnw:124-126
###################################################
data("npext", package="urca")   # load data
library("CADFtest")


###################################################
### code chunk number 3: CADFtest.Rnw:133-134
###################################################
ADFt <- CADFtest(npext$gnpperca, max.lag.y=3)


###################################################
### code chunk number 4: CADFtest.Rnw:139-140
###################################################
ADFt$p.value


###################################################
### code chunk number 5: CADFtest.Rnw:145-146
###################################################
CADFpvalues(ADFt$statistic, type="trend", rho2=1)


###################################################
### code chunk number 6: CADFtest.Rnw:153-154
###################################################
print(ADFt)


###################################################
### code chunk number 7: CADFtest.Rnw:161-162
###################################################
summary(ADFt)


###################################################
### code chunk number 8: CADFtest.Rnw:169-170
###################################################
res.ADFt <- residuals(ADFt)


###################################################
### code chunk number 9: CADFtest.Rnw:175-176
###################################################
plot(ADFt)


###################################################
### code chunk number 10: CADFtest.Rnw:181-182
###################################################
plot(ADFt)


###################################################
### code chunk number 11: CADFtest.Rnw:190-191
###################################################
plot(ADFt, plots=c(1,3,4))


###################################################
### code chunk number 12: CADFtest.Rnw:200-201
###################################################
plot(ADFt, plots=c(1,3,4))


###################################################
### code chunk number 13: CADFtest.Rnw:208-212
###################################################
npext$unemrate <- exp(npext$unemploy)      # compute unemployment rate
L <- ts(npext, start=1860)                 # time series of levels
D <- diff(L)                               # time series of diffs
S <- window(ts.intersect(L,D), start=1909) # select same sample as Hansen's


###################################################
### code chunk number 14: CADFtest.Rnw:219-220
###################################################
ADFt <- CADFtest(L.gnpperca ~ 1, data = S, max.lag.y = 3)


###################################################
### code chunk number 15: CADFtest.Rnw:227-228
###################################################
CADFtest(L.gnpperca ~ 1, data = S, max.lag.y = 4, criterion = "BIC", dname = "Extended Nelson-Plosser data")


###################################################
### code chunk number 16: CADFtest.Rnw:233-234
###################################################
ADFt2 <- update(ADFt, change=list("max.lag.y = 4", "criterion = 'BIC'", "dname = 'Extended Nelson-Plosser data'"))


###################################################
### code chunk number 17: CADFtest.Rnw:241-242
###################################################
ADFt <- CADFtest(L.gnpperca ~ 1, data = S, max.lag.y = 3)


###################################################
### code chunk number 18: CADFtest.Rnw:247-249
###################################################
CADFt <- update(ADFt, change=list("+ D.unemrate", "kernel = 'Parzen'", "prewhite = FALSE"))
print(CADFt)


###################################################
### code chunk number 19: CADFtest.Rnw:256-258
###################################################
CADFt <- update(CADFt, change=list("max.lag.X = 3", "min.lag.X = -3", "criterion = 'BIC'"))
print(CADFt)


###################################################
### code chunk number 20: CADFtest.Rnw:263-264
###################################################
CADFt <- CADFtest(L.gnpperca ~ D.unemrate, data = S, max.lag.y = 3, max.lag.X = 3, min.lag.X = -3, criterion = "BIC", kernel = "Parzen", prewhite = FALSE)


###################################################
### code chunk number 21: CADFtest.Rnw:269-270
###################################################
CADFt <- CADFtest(L.gnpperca ~ D.unemrate + D.indprod, data = S, max.lag.y = 3, max.lag.X = 3, min.lag.X = -3, criterion = "BIC", kernel = "Parzen", prewhite = FALSE)


###################################################
### code chunk number 22: CADFtest.Rnw:331-333
###################################################
CADFpvalues(t0 = -2.2, rho2 = 0.53)
CADFpvalues(t0 = -1.7, rho2 = 0.20)


###################################################
### code chunk number 23: CADFtest.Rnw:340-341
###################################################
CADFpvalues(-0.44, type = "drift", rho2 = 1)


###################################################
### code chunk number 24: CADFtest.Rnw:357-359
###################################################
library("urca")
adf.urca <- ur.df(npext$gnpperca[-(1:49)], type = "trend", lags = 3)


###################################################
### code chunk number 25: CADFtest.Rnw:366-368
###################################################
library("tseries")
adf.tseries <- adf.test(npext$gnpperca[-(1:49)], k = 3)


###################################################
### code chunk number 26: CADFtest.Rnw:375-377
###################################################
library("fUnitRoots")
adf.fUnitRoots <- unitrootTest(npext$gnpperca, lags = 3, type = "ct")



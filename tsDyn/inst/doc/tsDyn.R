### R code from vignette source 'tsDyn.Stex'
### Encoding: UTF-8

###################################################
### code chunk number 1: tsDyn.Stex:12-15
###################################################
require(tsDyn)
require(tseriesChaos)
options(prompt=" ", encoding="LATIN-9")


###################################################
### code chunk number 2: tsDyn.Stex:59-60 (eval = FALSE)
###################################################
## autopairs(x, lag=, type=, h=)


###################################################
### code chunk number 3: tsDyn.Stex:101-103
###################################################
obj <- llar(log(lynx), m=3)
plot(obj)


###################################################
### code chunk number 4: tsDyn.Stex:107-108
###################################################
obj <- data.frame(obj)


###################################################
### code chunk number 5: tsDyn.Stex:111-112
###################################################
names(obj)


###################################################
### code chunk number 6: tsDyn.Stex:115-116
###################################################
plot(RMSE~eps, data=obj, type="l", log="x")


###################################################
### code chunk number 7: tsDyn.Stex:128-129 (eval = FALSE)
###################################################
## delta.test(x)


###################################################
### code chunk number 8: tsDyn.Stex:153-154
###################################################
availableModels()


###################################################
### code chunk number 9: tsDyn.Stex:168-169 (eval = FALSE)
###################################################
## x.new <- predict(obj, n.ahead = )


###################################################
### code chunk number 10: tsDyn.Stex:206-208
###################################################
usual <- linear(lynx, m=3)
adf<- linear(lynx, m=2, type="ADF")


###################################################
### code chunk number 11: tsDyn.Stex:213-215
###################################################
 all.equal(deviance(adf), deviance(usual))
all.equal(residuals(usual), residuals(adf))


###################################################
### code chunk number 12: tsDyn.Stex:235-236 (eval = FALSE)
###################################################
## obj <- setar(x, m=, d=, steps=, thDelay= )


###################################################
### code chunk number 13: tsDyn.Stex:244-245 (eval = FALSE)
###################################################
## obj <- setar(x, m=, d=, steps=, mTh= )


###################################################
### code chunk number 14: tsDyn.Stex:249-250 (eval = FALSE)
###################################################
## obj <- setar(x, m=, d=, steps=, thVar= )


###################################################
### code chunk number 15: tsDyn.Stex:256-257 (eval = FALSE)
###################################################
## obj <- setar(x, m=, d=, steps=, thDelay=, nthresh=2)


###################################################
### code chunk number 16: tsDyn.Stex:269-270 (eval = FALSE)
###################################################
## obj <- setar(x, m=, d=, steps=, thDelay = , mL =, mH =)


###################################################
### code chunk number 17: tsDyn.Stex:274-275 (eval = FALSE)
###################################################
## obj <- setar(x, m=, d=, steps=, thDelay = , ML =, MH =)


###################################################
### code chunk number 18: tsDyn.Stex:310-311 (eval = FALSE)
###################################################
## obj <- nnetTs(x, m=, d=, steps=, size=)


###################################################
### code chunk number 19: tsDyn.Stex:324-325 (eval = FALSE)
###################################################
## obj <- aar(x, m=, d=, steps=)


###################################################
### code chunk number 20: tsDyn.Stex:334-336
###################################################
x <- log10(lynx)
selectSETAR(x, m=3, mL=1:3, mH=1:3, thSteps = 5, thDelay=0:2)


###################################################
### code chunk number 21: tsDyn.Stex:347-350
###################################################
str(lynx)
summary(lynx)
plot(lynx)


###################################################
### code chunk number 22: tsDyn.Stex:357-358
###################################################
x <- log10(lynx)


###################################################
### code chunk number 23: tsDyn.Stex:362-367
###################################################
par(mfrow=c(2,1), mar=c(0,0,0,0))
plot(x, ax=F)
box()
plot(x[length(x):1], type="l", ax=F)
box()


###################################################
### code chunk number 24: tsDyn.Stex:371-374
###################################################
par(mfrow=c(2,1), mar=c(2,2,0,0))
autopairs(x, lag=1, type="regression")
autopairs(x, lag=3, type="regression")


###################################################
### code chunk number 25: tsDyn.Stex:380-381
###################################################
hist(x, br=13)


###################################################
### code chunk number 26: tsDyn.Stex:385-388
###################################################
par(mfrow=c(2,1), mar=c(2,4,0,0))
acf(x)
pacf(x)


###################################################
### code chunk number 27: tsDyn.Stex:392-394
###################################################
library(tseriesChaos)
mutual(x)


###################################################
### code chunk number 28: tsDyn.Stex:398-399
###################################################
recurr(x, m=3, d=1, levels=c(0,0.2,1))


###################################################
### code chunk number 29: tsDyn.Stex:405-406
###################################################
lag.plot(x, lags=3, layout=c(1,3))


###################################################
### code chunk number 30: tsDyn.Stex:412-414
###################################################
delta.test(x)
delta.lin.test(x)


###################################################
### code chunk number 31: tsDyn.Stex:428-430
###################################################
mod.ar <- linear(x, m=2)
mod.ar


###################################################
### code chunk number 32: tsDyn.Stex:434-436
###################################################
mod.setar <- setar(x, m=2, mL=2, mH=2, thDelay=1)
mod.setar


###################################################
### code chunk number 33: tsDyn.Stex:439-440
###################################################
beta <- round(coef(mod.setar),3)


###################################################
### code chunk number 34: tsDyn.Stex:451-452
###################################################
set.seed(10)


###################################################
### code chunk number 35: tsDyn.Stex:454-460
###################################################
mod <- list()
mod[["linear"]] <- linear(x, m=2)
mod[["setar"]] <- setar(x, m=2, thDelay=1)
mod[["lstar"]] <- lstar(x, m=2, thDelay=1)
mod[["nnetTs"]] <- nnetTs(x, m=2, size=3)
mod[["aar"]] <- aar(x, m=2)


###################################################
### code chunk number 36: tsDyn.Stex:465-467
###################################################
sapply(mod, AIC)
sapply(mod, MAPE)


###################################################
### code chunk number 37: tsDyn.Stex:473-474
###################################################
summary(mod[["setar"]])


###################################################
### code chunk number 38: tsDyn.Stex:478-479 (eval = FALSE)
###################################################
## plot(mod[["setar"]])


###################################################
### code chunk number 39: tsDyn.Stex:485-494
###################################################
set.seed(10)
mod.test <- list()
x.train <- window(x, end=1924)
x.test <- window(x, start=1925)
mod.test[["linear"]] <- linear(x.train, m=2)
mod.test[["setar"]] <- setar(x.train, m=2, thDelay=1)
mod.test[["lstar"]] <- lstar(x.train, m=2, thDelay=1, trace=FALSE, control=list(maxit=1e5))
mod.test[["nnet"]] <- nnetTs(x.train, m=2, size=3, control=list(maxit=1e5))
mod.test[["aar"]] <- aar(x.train, m=2)


###################################################
### code chunk number 40: tsDyn.Stex:498-504
###################################################
frc.test <- lapply(mod.test, predict, n.ahead=10)
plot(x.test,ylim=range(x))
for(i in 1:length(frc.test))
	lines(frc.test[[i]], lty=i+1, col=i+1)

legend(1925,2.4, lty=1:(length(frc.test)+1), col=1:(length(frc.test)+1), legend=c("observed",names(frc.test)))


###################################################
### code chunk number 41: tsDyn.Stex:510-511
###################################################
par(cex=0.6)


###################################################
### code chunk number 42: tsDyn.Stex:516-518
###################################################
x.new <- predict(mod[["linear"]], n.ahead=100)
lag.plot(x.new, 1)


###################################################
### code chunk number 43: tsDyn.Stex:522-524
###################################################
x.new <- predict(mod[["setar"]], n.ahead=100)
lag.plot(x.new, 1)


###################################################
### code chunk number 44: tsDyn.Stex:529-531
###################################################
x.new <- predict(mod[["nnetTs"]], n.ahead=100)
lag.plot(x.new, 1)


###################################################
### code chunk number 45: tsDyn.Stex:544-546
###################################################
mod.point <- setar(x, m=10, mL=3, mH=10, thDelay=0, th=3.12)
lag.plot(predict(mod.point, n.ahead=100))


###################################################
### code chunk number 46: tsDyn.Stex:550-552
###################################################
mod.unstable <- setar(x, m=9, mL=9, mH=6, thDelay=4, th=2.61)
lag.plot(predict(mod.unstable, n.ahead=100))


###################################################
### code chunk number 47: tsDyn.Stex:556-558
###################################################
mod.chaos1 <- setar(x, m=5, mL=5, mH=3, thDelay=1, th=2.78)
lag.plot(predict(mod.chaos1, n.ahead=100))


###################################################
### code chunk number 48: tsDyn.Stex:560-562
###################################################
mod.chaos2 <- setar(x, m=5, mL=5, mH=3, thDelay=1, th=2.95)
lag.plot(predict(mod.chaos2, n.ahead=100))


###################################################
### code chunk number 49: tsDyn.Stex:576-581
###################################################
N <- 1000
x.new <- predict(mod[["setar"]], n.ahead=N)
x.new <- x.new + rnorm(N, sd=sd(x.new)/100)
ly <- lyap_k(x.new, m=2, d=1, t=1, k=2, ref=750, s=200, eps=sd(x.new)/10)
plot(ly)


###################################################
### code chunk number 50: tsDyn.Stex:585-589
###################################################
x.new <- predict(mod.chaos2, n.ahead=N)
x.new <- x.new + rnorm(N, sd=sd(x.new)/100)
ly <- lyap_k(x.new, m=5, d=1, t=1, k=2, ref=750, s=200, eps=sd(x.new)/10)
plot(ly)


###################################################
### code chunk number 51: tsDyn.Stex:593-594
###################################################
lyap(ly,start=6,end=70)


###################################################
### code chunk number 52: tsDyn.Stex:711-714 (eval = FALSE)
###################################################
## plot(setar(lynx, m=1, model="MTAR"))
## plot(setar(lynx, m=3, type="ADF"))
##  #both won't work



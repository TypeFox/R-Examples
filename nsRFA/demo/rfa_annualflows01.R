# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #
#       REGIONAL FREQUENCY ANALYSIS OF ANNUAL FLOWS IN PIEMONTE AND VALLE D'AOSTA       #
# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #

# ------------------------------------------------------------------------------------- #
#                         Regionalization of the index-flow                             #
# ------------------------------------------------------------------------------------- #

#readline()


# Data loading

data(hydroSIMN)

str(parameters)


#readline()

# Catchment parameters plot

plot(parameters)
par(ask = interactive())

#readline()

# Catchment parameters plot

op <- par(mfrow=c(2,2))
plot(parameters[c("Xbar","Ybar")])
plot(parameters[c("Hm","Ybar")])
plot(parameters[c("Xbar","S")])
plot(parameters[c("Hm","S")])
par(op)


#readline()

# Choice of the best regressions:
#  create a function using the 'leaps' function of package 'subselect'
#  to perform all-possible-regressions
# bestregressions <- function(dip,ind) {
#  Y <- as.numeric(dip)
#  X <- ind
#  Sy <- var(Y)
#  Sx <- var(X)
#  Sxy <- var(X,Y)
#  Dm.mat <- Sx
#  Dm.H <- Sxy %*% t(Sxy)/Sy
#  require(subselect)
#  Dm.leaps <- leaps(Dm.mat, kmin=1, kmax=3, H=Dm.H, r=1, nsol=3)
#  Dm.leaps
#  for(i in 3:1) {for(j in 1:3) {print(colnames(X)[Dm.leaps$subsets[j,c(1:3),i]])}}
# }

dataregr <- cbind(parameters,(parameters["Dm"])^(1/3),log(parameters["Dm"]),log(parameters["Am"]))
names(dataregr) <- c(names(parameters),"sq3.Dm","ln.Dm","ln.Am")

# Regressions between Dm and parameters
# bestregressions(dataregr[,2],dataregr[,-c(1,2,17,18)])
# [1] "Am"    "S2000" "IT"
# [1] "S2000" "IT"    "ln.Am"
# [1] "Am"    "S2000" "IB"
# [1] "S2000" "ln.Am"
# [1] "Am"    "S2000"
# [1] "Hm"    "ln.Am"
# [1] "IT"
# [1] "IB"
# [1] "ln.Am"


#readline()

# Regressions between sq3.Dm and parameters
# bestregressions(dataregr[,17],dataregr[,-c(1,2,17,18)])
# [1] "Hm"   "NORD" "IB"
# [1] "S2000" "IT"    "ln.Am"
# [1] "Hm"    "S2000" "ln.Am"
# [1] "S2000" "ln.Am"
# [1] "Hm"    "ln.Am"
# [1] "Hm" "IB"
# [1] "IT"
# [1] "IB"
# [1] "ln.Am"


#readline()

# Regressions between ln.Dm and parameters
# bestregressions(dataregr[,18],dataregr[,-c(1,2,17,18)])
# [1] "Hm"   "NORD" "IB"
# [1] "Hm"    "NORD"  "ln.Am"
# [1] "Hm"    "EST"   "ln.Am"
# [1] "S2000" "ln.Am"
# [1] "Hm"    "ln.Am"
# [1] "Hm" "IB"
# [1] "IT"
# [1] "IB"
# [1] "ln.Am"


#readline()

# Best regression

regr01 <- lm(ln.Dm ~ 1 + Hm + NORD + IB, dataregr)
print(summary(regr01))

print(R2.lm(regr01))

print(prt.lm(regr01))

print(vif.lm(regr01))

# RMSE
resid <- exp(regr01$fitted.values) - parameters[,c("Dm")]
print(sqrt(sum((resid)^2)/length(resid)))

# RMSEcv
predictions <- jackknife1.lm(regr01)
resid <- exp(predictions) - parameters[,c("Dm")]
print(sqrt(sum((resid)^2)/length(resid)))

print(round(cor(regr01$model[-1]),3))


#readline()

# Best regression plots

op <- par(mfrow=c(2,2))
plot(regr01$fitted.values,regr01$residuals,xlab="Fitted",ylab="Residuals")
abline(0,0,lty=3)
normplot(regr01$residuals,xlab="Residuals")
plot(parameters[,c("Dm")],exp(regr01$fitted.values),xlab="Originals",ylab="Fitted")
abline(0,1,lty=3)
intervals <- predinterval.lm(regr01)
intervals <- intervals[order(intervals[,1]),]
plot(parameters[,c("Dm")],exp(predictions),xlab="Originals",ylab="Predicted")
abline(0,1,lty=3)
lines(exp(intervals[,c(1,2)]),lty=2)
lines(exp(intervals[,c(1,3)]),lty=2)
par(op)


#readline()

# Simplest regression

regr02 <- lm(sq3.Dm ~ 1 + ln.Am + Hm, dataregr)
print(summary(regr02))

print(R2.lm(regr02))

print(prt.lm(regr02))

print(vif.lm(regr02))

# RMSE
resid <- regr02$fitted.values^3 - parameters[,c("Dm")]
print(sqrt(sum((resid)^2)/length(resid)))

# RMSEcv
predictions <- jackknife1.lm(regr02)
resid <- predictions^3 - parameters[,c("Dm")]
print(sqrt(sum((resid)^2)/length(resid)))

print(round(cor(regr02$model[-1]),3))


#readline()

# Best regression plots

op <- par(mfrow=c(2,2))
plot(regr02$fitted.values,regr02$residuals,xlab="Fitted",ylab="Residuals")
abline(0,0,lty=3)
normplot(regr02$residuals,xlab="Residuals")
plot(parameters[,c("Dm")],(regr02$fitted.values)^3,xlab="Originals",ylab="Fitted")
abline(0,1,lty=3)
intervals <- predinterval.lm(regr02)
intervals <- intervals[order(intervals[,1]),]
plot(parameters[,c("Dm")],predictions^3,xlab="Originals",ylab="Predicted")
abline(0,1,lty=3)
lines((intervals[,c(1,2)])^3,lty=2)
lines((intervals[,c(1,3)])^3,lty=2)
par(op)


# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #
#                                        THE END                                        #
# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #




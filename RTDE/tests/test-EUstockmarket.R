data(EuStockMarkets)


library(RTDE)
library(tseries)

FTSE <- diff(log(EuStockMarkets))[,"FTSE"]
CAC <- diff(log(EuStockMarkets))[,"CAC"]
DAX <- diff(log(EuStockMarkets))[,"DAX"]

ftse.garch <- garch(FTSE)
cac.garch <- garch(CAC)
dax.garch <- garch(DAX)

X <- diff(as.numeric(dax.garch$residuals))
X <- X[!is.na(X)]
Y <- diff(as.numeric(cac.garch$residuals))
Y <- Y[!is.na(Y)]

stockreturn <- dataRTDE(cbind(X, Y))
stockZ <- zvalueRTDE(cbind(X,Y), omega=1/2, length(X)-1, output="relexcess")

plot(stockreturn, which=1)
plot(stockreturn, which=2)
plot(apply(stockreturn$data, 2, function(x) rank(x)/length(x)))
qqparetoplot(stockZ$Z)


feta <- RTDE(cbind(X, Y), nbpoint=seq(100, 300, by=100), alpha=c(0, 0.5), omega=1/2:3)
feta2 <- RTDE(cbind(X, Y), nbpoint=seq(100, 300, by=100), alpha=c(0, 0.5), omega=1/2)
summary(feta)
# plot(feta, which=1)

fprob <- prob(feta, q=5:6*100)
plot(fprob, which=3)

fprob <- prob(feta, q=500)

plot(fprob, which=3)

fprob <- prob(feta2, q=500)

plot(fprob, which=3)


feta <- RTDE(cbind(X, Y), nbpoint=seq(10, 240, by=10), alpha=c(0, 0.25, 0.5), omega=1/2)
fprob <- prob(feta, q=5)
plot(fprob, which=3)


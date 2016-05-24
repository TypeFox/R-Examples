mleurDiag <-
function(y, lag.max="default"){
n <- length(y)
if (identical(lag.max, "default") || lag.max < 2) lag.max <- max(min(n/4, 20), 2)
yc <- as.vector(y-mean(y))
phiH <- ar1est(yc) #default is mean correction now
res <- yc[2:n]-phiH*yc[1:(n-1)]
r <- (acf(res, lag.max = lag.max, plot=FALSE)$acf)[-1]
Q <- (n*(n+2))*cumsum((r^2)/(n-(1:lag.max)))
Q <- Q[-1]
pv <- 1-pchisq(Q, df=1:(lag.max-1))
A <- (phiH^2)^(1:(lag.max-1))
B <- A*phiH^2
moe <- 1.96*sqrt(c(phiH^2,1-A+B)/n)
rmax <- max(abs(r), moe)
par(mfrow=c(3, 1))
boxplot(res, horizontal=TRUE, notch=TRUE,cex=2)
out<-shapiroTest(res)
pvW<-as.numeric(out@test$p.value)
title(sub=paste("Wilk-Shapiro test, p-value: ", format.pval(pvW)))
plot(2:lag.max, pv, abline(h=0.05, col="red", lty=2), ylim=c(0,1),
    xlab="lag", ylab="p-value", pch=16)
title(main="P-values of Box-Ljung test")
plot(1:lag.max, r, type="h", lwd = 2, xlab="lag", ylab="acf", ylim=c(-rmax, rmax))
title(main="Residual autocorrelation plot")
abline(h=0, col="gray")
lines(1:lag.max, moe, col="red", lty=2)
lines(1:lag.max, -moe, col="red", lty=2)
par(mfrow=c(1,1))
invisible(res)
}

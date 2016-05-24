## ----message=FALSE-------------------------------------------------------
library(bcpa)

## ----GetRhoDemo1, out.width="\\textwidth", fig.height=3, echo=-1---------
par(bty="l")
rho <- 0.8
x.full <- arima.sim(1000, model=list(ar = rho))
t.full <- 1:1000
keep <- sort(sample(1:1000, 200))
x <- x.full[keep]
t <- t.full[keep]
plot(t,x, type="l")

## ----GetRhoDemo2, out.width="\\textwidth", fig.height=3, echo=-1, tidy=FALSE----
par(bty="l")
rhos <- seq(0,.99,.01)
L <- rep(NA, length(rhos))
for(i in 1:length(rhos))
  L[i] <- GetL(x,t,rhos[i])
# plot likelihood profile
plot(rhos, L, type="l")
abline(v = rho, lty=2, lwd=2); abline(v = rhos[L == max(L)], lty=3, lwd=2)
legend("bottomleft", legend=c("true value","MLE"), lty=2:3, lwd=2)

## ------------------------------------------------------------------------
GetRho(x, t, tau=FALSE)
GetRho(x, t, tau=TRUE)

## ----Simp, out.width="0.7\\textwidth", fig.width=5, fig.height=5, echo=-1----
par(bty="l", cex.lab=1.25)
data(Simp)
head(Simp)
plot(Simp)

## ----MakeTrack-----------------------------------------------------------
X <- cumsum(arima.sim(n=100, model=list(ar=0.8)))
Y <- cumsum(arima.sim(n=100, model=list(ar=0.8)))
Time <- 1:100
mytrack <- MakeTrack(X,Y,Time)
plot(mytrack)

## ------------------------------------------------------------------------
Simp.VT <- GetVT(Simp)
head(Simp.VT)

## ----Histograms, fig.height = 4, out.width="0.8\\textwidth"--------------
par(mfrow=c(1,2))
hist(Simp.VT$V, breaks=20, col="grey")
hist(Simp.VT$Theta, breaks=20, col="grey")

## ----OneBreak, fig.height=4, echo=-1, tidy=FALSE, cache=TRUE-------------
set.seed(2)
par(bty="l")
mu1 <- 5; mu2 <- 3
sigma1 <- 2; sigma2 <- 1
rho1 <- 0.5; rho2 <- 0.5

SimTS <- function(n, mu, rho, sigma)
{
  X.standard <- arima.sim(n, model=list(ar = rho))
  X.standard/sd(X.standard)*sigma + mu
}

# create time series with break at 500
t.full <- 1:1000
t.break <- 500
x.full <- c(SimTS(t.break, mu1, rho1, sigma1), 
            SimTS(max(t.full)-t.break+1, mu2, rho2, sigma2))

# subsample 100 observations and estimate
keep <- sort(sample(1:length(x.full), 100))
x <- x.full[keep]
t <- t.full[keep]
(BB <- GetBestBreak(x,t, tau=FALSE))

## ----OneBreak2, fig.height=4, echo=-1, tidy=FALSE, out.width="0.8\\textwidth"----
par(bty="l")
plot(t,x, type="l")
abline(v = 500, col=2, lwd=2, lty=2); abline(v = BB[2], col=2, lwd=2, lty=3)
legend("topright", legend=c("true break", "estimated break"), col=2, lwd=2, lty=2:3)

## ------------------------------------------------------------------------
GetModels(x,t,BB[1], tau=FALSE)

## ------------------------------------------------------------------------
GetModels(x,t,BB[1], tau=FALSE, K=0.5)

## ------------------------------------------------------------------------
GetModels(x,t,BB[1], tau=FALSE, K=5)

## ----OneBreak3, fig.height=4, echo=1:3, out.width="0.8\\textwidth"-------
mu1 <- 0; mu2 <- 0
sigma1 <- 1; sigma2 <- 1
rho1 <- 0.9; rho2 <- 0.2
set.seed(11)
SimTS <- function(n, mu, rho, sigma)
{
  X.standard <- arima.sim(n, model=list(ar = rho))
  X.standard/sd(X.standard)*sigma + mu
}

# create time series with break at 500
t.full <- 1:1000
t.break <- 500
x.full <- c(SimTS(t.break, mu1, rho1, sigma1), SimTS(max(t.full)-t.break+1, mu2, rho2, sigma2))

# subsample 100 observations and estimate
keep <- sort(sample(1:length(x.full), 100))
x <- x.full[keep]
t <- t.full[keep]
BB <- GetBestBreak(x,t, tau=FALSE)

par(bty="l")
plot(t,x, type="l")

abline(v = 500, col=2, lwd=2, lty=2)
abline(v = BB[2], col=2, lwd=2, lty=3)
legend("topright", legend=c("true break", "estimated break"), col=2, lwd=2, lty=2:3)

## ------------------------------------------------------------------------
GetModels(x,t,BB[1], tau=FALSE)

## ----LampreyBCPA, cache=TRUE---------------------------------------------
  Simp.ws <- WindowSweep(Simp.VT, "V*cos(Theta)", windowsize=50, progress=FALSE, K=2)

## ------------------------------------------------------------------------
head(Simp.ws$ws) 

## ----BCPAsmooth, fig.width = 9, fig.height=5, out.width="\\textwidth", size="small", echo=-1, warning=FALSE----
par(mfrow=c(2,1), mar=c(0,4,0,1), oma=c(4,0,1,0), bty="l")
plot(Simp.ws, type="smooth")
plot(Simp.ws, type="smooth", threshold = 7)

## ----BCPAflat, fig.width = 9, fig.height=5, out.width="\\textwidth", size="small", echo=-1, warning=FALSE----
par(mfrow=c(2,1), mar=c(0,4,0,1), oma=c(4,0,1,0), bty="l")
plot(Simp.ws, type="flat")
plot(Simp.ws, type="flat", clusterwidth=3)

## ------------------------------------------------------------------------
ChangePointSummary(Simp.ws, clusterwidth=3)

## ----BCPApaths, echo=-1, warning=FALSE, fig.height=4, fig.width=8--------
par(mfrow=c(1,2))
PathPlot(Simp, Simp.ws, type="flat", clusterwidth = 3, main="Flat BCPA")
PathPlot(Simp, Simp.ws, type="smooth", main="Smooth BCPA")

## ----out.width="0.6\\textwidth", fig.height=4, fig.width=5, echo=-1------
par(bty="l")
PhasePlot(Simp.ws, type="smooth", clusterwidth = 3)

## ----DiagnosticPlot, fig.width = 8, fig.height=3, out.width = "0.8\\textwidth", echo=-1----
par(bty="l")
DiagPlot(Simp.ws)


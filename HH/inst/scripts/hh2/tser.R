### R code from vignette source '~/WindowsC/HOME/rmh/hh.e2/hh2/tser.tex'

###################################################
### code chunk number 1: tser.tex:9-10
###################################################
library(HH)


###################################################
### code chunk number 2: tser.tex:118-128
###################################################
## hhcapture("tseralign1.Rout", '
x <- ts(sample(10), start=1978)
y <- ts(sample(6), start=1980)
x
y
ts.union(x,y)
## ')
## hhcapture("tseralign2.Rout", '
ts.intersect(x,y)
## ')


###################################################
### code chunk number 3: tser.tex:390-415
###################################################
## hhcapture("my.pacf.Rout", '
## This function illustrates the definition of the pacf.
## Do NOT use in actual calculations!
##
## my.pacf requires a detrended series, otherwise the answer is
## nonsense, as it starts losing precision after the first few lags.

my.pacf <- function(z, k=2) {
  z <- z - mean(z)
  x <- ts.intersect(z, lag(z,-1))
  if (k==1) return(cor(x[,1], x[,2]))
  for (kk in 2:k) x <- ts.intersect(x, lag(z,-kk))
  nr <- nrow(x)
  nc <- ncol(x)
  r1 <- lm(x[,1]  ~ -1 + x[,-c(1,nc)])$resid
  r2 <- lm(x[,nc] ~ -1 + x[,-c(1,nc)])$resid
  cor(r1,r2)
}

data(ozone)
## 8/31 -- 11/53 has no missing data.
ozone.subset <- ts(ozone[68:335], start=c(1931,8), freq=12)
my.pacf(ozone.subset, 2)
acf(ozone.subset, type="partial", plot=FALSE)$acf[2]
## ')


###################################################
### code chunk number 4: tser.tex:800-806
###################################################
data(elnino)
## hhpdf("elnino1.pdf", height=8, width=7)
print(portrait=TRUE,
tsacfplots(elnino, main="", cex=1.3, ylab="elnino --- Southern Oscillation", lwd=1)
)
## hhdev.off()


###################################################
### code chunk number 5: tser.tex:819-824
###################################################
## hhpdf("elnino2.pdf", height=8, width=7)
print(portrait=TRUE,
tsacfplots(diff(elnino), main="", cex=1.3, lwd=1)
)
## hhdev.off()


###################################################
### code chunk number 6: tser.tex:859-871
###################################################
## hhcapture("elnino3.Rout", '
elnino.loop <- arma.loop(elnino, order=c(2,1,2))
elnino.loop
## ')
## hhpdf("elnino3.pdf", height=13, width=9)
elnino.diags <- diag.arma.loop(elnino.loop, lag.max=24)
print(portrait=TRUE,
tsdiagplot(armas=elnino.loop, diags=elnino.diags,
	   lag.lim=c(-2,22),
           main="", lwd=1)
)
## hhdev.off()


###################################################
### code chunk number 7: tser.tex:1099-1104
###################################################
## hhpdf("tsamsta1.pdf", height=8, width=7)
print(portrait=TRUE,
tsacfplots(co2,                   xlab=NULL, ylab=list("co2", cex=1.5),                   lag.max=36, cex=1.3, main="", lwd=1)
)
## hhdev.off()


###################################################
### code chunk number 8: tser.tex:1119-1123
###################################################
## hhpdf("co2_10year.pdf", height=4, width=7)
seqplot(ts(co2[1:120], start=c(1959,1), frequency=12),
        xlab=NULL, ylab=list("co2", cex=1.5), cex=1.5, main="", lwd=1)
## hhdev.off()


###################################################
### code chunk number 9: tser.tex:1148-1153
###################################################
## hhpdf("tsamsta1a.pdf", height=8, width=7)
print(portrait=TRUE,
tsacfplots(diff(co2,1),           xlab=NULL, ylab=list("diff(co2,1)", cex=1.5),           lag.max=36, cex=1.3, main="", lwd=1)
)
## hhdev.off()


###################################################
### code chunk number 10: tser.tex:1174-1179
###################################################
## hhpdf("tsamsta2.pdf", height=8, width=7)
print(portrait=TRUE,
tsacfplots(diff(diff(co2,1), 12), xlab=NULL, ylab=list("diff(diff(co2,1), 12)", cex=1.5), lag.max=36, cex=1.3, main="", lwd=1)   ## black and white
)
## hhdev.off()


###################################################
### code chunk number 11: tser.tex:1191-1196
###################################################
## not included in printed book
##HT main="Figure 2alt (not in printed paper).  diff(diff(co2,1), 12)  with lag.0=FALSE"
print(portrait=TRUE,
tsacfplots(diff(diff(co2,1), 12), xlab=NULL, ylab=list("diff(diff(co2,1), 12)", cex=1.5), lag.max=36, cex=1.3, main="", lwd=1, lag.0=FALSE)
)


###################################################
### code chunk number 12: tser.tex:1235-1250
###################################################
ddco2.loop <- arma.loop(co2,
                        order=c(2,1,2),
                        seasonal=list(order=c(0,1,0), period=12))
ddco2.loop

ddco2.diags <- diag.arma.loop(ddco2.loop, co2, lag.max=36)
## hhpdf("tsamsta3.pdf", height=13, width=9)
print(portrait=TRUE,
tsdiagplot(armas=ddco2.loop, diags=ddco2.diags,
	   lag.lim=c(-2,38),
	   lag.x.at=seq(0,36,6),
	   lag.x.labels=c(0,"",12,"",24,"",36),
           main="", lwd=1)
)
## hhdev.off()


###################################################
### code chunk number 13: tser.tex:1310-1319
###################################################
old.width <- options(width=70)
## hhcapture("tsamsta4.Rout", '
ddco2.loopPQ <-
       arma.loop(co2,
                 order=c(2,1,2),
                 seasonal=list(order=c(0,1,1), period=12))
ddco2.loopPQ
## ')
options(old.width)


###################################################
### code chunk number 14: tser.tex:1332-1342
###################################################
ddco2.diagsPQ <- diag.arma.loop(ddco2.loopPQ, co2, lag.max=36)
## hhpdf("tsamsta4.pdf", height=13, width=9)
print(portrait=TRUE,
tsdiagplot(armas=ddco2.loopPQ, diags=ddco2.diagsPQ,
	   lag.lim=c(-2,38),
	   lag.x.at=seq(0,36,6),
	   lag.x.labels=c(0,"",12,"",24,"",36),
           main="", lwd=1)
)
## hhdev.off()


###################################################
### code chunk number 15: tser.tex:1427-1434
###################################################
## hhcapture("tsamsta011011.Rout", '
co2.arima <- ddco2.loopPQ[["0","1"]]
co2.coef.t <- co2.arima$coef / sqrt(diag(co2.arima$var.coef))
co2.arima
co2.coef.t
vcov(co2.arima)
## ')


###################################################
### code chunk number 16: tser.tex:1447-1454
###################################################
## hhcapture("tsamsta111011.Rout", '
co2.arima11 <- ddco2.loopPQ[["1","1"]]
co2.coef11.t <- co2.arima11$coef / sqrt(diag(co2.arima11$var.coef))
co2.arima11
co2.coef11.t
vcov(co2.arima11)
## ')


###################################################
### code chunk number 17: tser.tex:1476-1490
###################################################
## hhpdf("tsamsta5.pdf", height=3, width=7)
co2.forecast <- predict(co2.arima11, n.ahead=12)

co2.last.year <- ts(co2[457:468],
                    start=time(co2)[457],
                    frequency=frequency(co2))

seqplotForecast(co2.last.year, co2.forecast,
                x.at=seq(1997,1999,.5),
                x.labels=c(1997,"",1998,"",1999),
                xlim=c(1996.9,1999),
                ylim=c(359, 370),
                main="", lwd=1, cex=2, lty=c(1,2))
## hhdev.off()


###################################################
### code chunk number 18: tser.tex:1605-1630
###################################################
## Mystery time series X
data(tser.mystery.X)

## hhpdf("arima_sim_X.pdf", height=8, width=7)
print(portrait=TRUE,
tsacfplots(tser.mystery.X, lwd=1, pch.seq=16, cex=.7)
)
## hhdev.off()

X.loop <- arma.loop(tser.mystery.X, order=c(2,0,2))
X.diag <- rearrange.diag.arma.loop(diag.arma.loop(X.loop, x=tser.mystery.X))

## hhpdf("arima_sim_Xd.pdf", height=13, width=9)
print(portrait=TRUE,
tsdiagplot(armas=X.loop, ts.diag=X.diag, lwd=1)
)
## hhdev.off()

## hhcapture("arima_sim_X.Rout", '
X.loop
## ')

## hhcapture("arima_sim_X11.Rout", '
X.loop[["1","1"]]
## ')


###################################################
### code chunk number 19: tser.tex:1695-1720
###################################################
## Mystery time series Y
data(tser.mystery.Y)

## hhpdf("arima_sim_Y.pdf", height=8, width=7)
print(portrait=TRUE,
tsacfplots(tser.mystery.Y, lwd=1, pch.seq=16, cex=.7)
)
## hhdev.off()

Y.loop <- arma.loop(tser.mystery.Y, order=c(2,0,2))
Y.diag <- rearrange.diag.arma.loop(diag.arma.loop(Y.loop, x=tser.mystery.Y))

## hhpdf("arima_sim_Yd.pdf", height=13, width=9)
print(portrait=TRUE,
tsdiagplot(armas=Y.loop, ts.diag=Y.diag, lwd=1)
)
## hhdev.off()

## hhcapture("arima_sim_Y.Rout", '
Y.loop
## ')

## hhcapture("arima_sim_Y11.Rout", '
Y.loop[["1","1"]]
## ')


###################################################
### code chunk number 20: tser.tex:1788-1813
###################################################
## Mystery time series Z
data(tser.mystery.Z)

## hhpdf("arima_sim_Z.pdf", height=8, width=7)
print(portrait=TRUE,
tsacfplots(tser.mystery.Z, lwd=1, pch.seq=16, cex=.7)
)
## hhdev.off()

Z.loop <- arma.loop(tser.mystery.Z, order=c(2,0,2))
Z.diag <- rearrange.diag.arma.loop(diag.arma.loop(Z.loop, x=tser.mystery.Z))

## hhpdf("arima_sim_Zd.pdf", height=13, width=9)
print(portrait=TRUE,
tsdiagplot(armas=Z.loop, ts.diag=Z.diag, lwd=1)
)
## hhdev.off()

## hhcapture("arima_sim_Z.Rout", '
Z.loop
## ')

## hhcapture("arima_sim_Z11.Rout", '
Z.loop[["1","1"]]
## ')


###################################################
### code chunk number 21: tser.tex:1938-1944
###################################################
data(product)
## hhpdf("prodfig1.pdf", height=8, width=7)
print(portrait=TRUE,
tsacfplots(product, main="", lwd=1)
)
## hhdev.off()


###################################################
### code chunk number 22: tser.tex:1958-1963
###################################################
## hhpdf("prodfig2.pdf", height=8, width=7)
print(portrait=TRUE,
tsacfplots(diff(product), main="", lwd=1)
)
## hhdev.off()


###################################################
### code chunk number 23: tser.tex:1985-1990
###################################################
## hhcapture("product.Rout", '
product.loop <- arma.loop(product, order=c(2,1,2))
product.diags <- diag.arma.loop(product.loop, x=product, lag.max=60)
product.loop
## ')


###################################################
### code chunk number 24: tser.tex:2002-2010
###################################################
## hhpdf("prodfig4.pdf", height=13, width=9)
print(portrait=TRUE,
tsdiagplot(armas=product.loop, diags=product.diags,
           x=product,
	   lag.lim=c(-2,60),
	   lag.x.at=seq(0,52,13), lwd=1)
)
## hhdev.off()


###################################################
### code chunk number 25: tser.tex:2056-2064
###################################################
## hhpdf("nottema.pdf", height=8, width=7)
print(portrait=TRUE,
tsacfplots(nottem, lwd=1, cex=1.3)
)
## hhdev.off()
print(portrait=TRUE,
tsacfplots(diff(nottem), lwd=1, cex=1.3)  ## not displayed in book
)


###################################################
### code chunk number 26: tser.tex:2074-2082
###################################################
## hhpdf("nottemb.pdf", height=8, width=7)
print(portrait=TRUE,
tsacfplots(diff(nottem, 12), lwd=1, cex=1.3)
)
## hhdev.off()
print(portrait=TRUE,
tsacfplots(diff(diff(nottem), 12), lwd=1, cex=1.3)  ## not displayed in book
)


###################################################
### code chunk number 27: tser.tex:2095-2105
###################################################
old.width <- options(width=70)
## hhcapture("nottem-3x3.Rout", '
nottem.loop <- arma.loop(nottem,  order=c(2,0,2),
                         seasonal=list(order=c(2,1,0), period=12),
                         method="ML")
print(nottem.loop, digits=4)
nottem.diag <-
   rearrange.diag.arma.loop(diag.arma.loop(nottem.loop, nottem))
## ')
options(old.width)


###################################################
### code chunk number 28: tser.tex:2117-2122
###################################################
## hhpdf("nottemc.pdf", height=13, width=9)
print(portrait=TRUE,
tsdiagplot(armas=nottem.loop, ts.diag=nottem.diag, lwd=1)
)
## hhdev.off()


###################################################
### code chunk number 29: tser.tex:2133-2136
###################################################
## hhcapture("nottem-100x210.Rout", '
nottem.loop[["1","0"]]
## ')


###################################################
### code chunk number 30: tser.tex:2192-2198
###################################################
data(tsq)
## hhpdf("tsq1.pdf", height=8, width=7)
print(portrait=TRUE,
tsacfplots(tsq, main="tsq 1", lwd=1)
)
## hhdev.off()


###################################################
### code chunk number 31: tser.tex:2220-2224
###################################################
## hhcapture("tsq.Rout", '
tsq.loop <- arma.loop(tsq, order=c(2,0,2))
tsq.loop
## ')


###################################################
### code chunk number 32: tser.tex:2249-2256
###################################################
## hhpdf("tsq2.pdf", height=13, width=9)
tsq.diag <-
  rearrange.diag.arma.loop(diag.arma.loop(tsq.loop, tsq))
print(portrait=TRUE,
tsdiagplot(armas=tsq.loop, ts.diag=tsq.diag, lwd=1)
)
## hhdev.off()


###################################################
### code chunk number 33: tser.tex:2295-2304
###################################################
data(ozone)
seqplot(ozone, "Missing values prevent calculation of acf and pacf.", cex=1.3) ## not in book
## 8/31 -- 11/53 has no missing data.
ozone.subset <- ts(ozone[68:335], start=c(1931,8), freq=12)
## hhpdf("ex0304-4.pdf", height=8, width=7)
print(portrait=TRUE,
tsacfplots(ozone.subset, cex=1.3, lwd=1)
)
## hhdev.off()


###################################################
### code chunk number 34: tser.tex:2337-2354
###################################################
acf.1  <- c(.8, .61, .47,  .40, .31,  .21, .18, .11, .06, .01)
pacf.1 <- c(.8, .08, .00, -.11, .00, -.12, .07, .05, .01, .02)
tmp <- data.frame(acf=c(1,acf.1, 1,pacf.1),
                  lag=c(0:10, 0:10),
                  type=rep(c("acf","pacf"),c(11,11)))
## hhpdf("co697ts1.pdf", width=7, height=3.5)
xyplot(acf ~ lag | type, data=tmp,
       n.used=100,
       panel=panel.acf, lwd=2,
       par.strip.text=list(cex=1.5),
       strip=function(...) strip.default(..., style=1),
       between=list(x=1, y=1),
       scales=list(alternating=FALSE, x=list(cex=1.2), y=list(cex=1.2)),
       xlab=list(cex=1.2), ylab=list("",cex=1.2),
       main=list("Time Series Question, n=100, Z.bar=25", cex=1.6),
       layout=c(2,1))
## hhdev.off()


###################################################
### code chunk number 35: tser.tex:2388-2406
###################################################
z22  <- c( .93, .92,.90,.90, .87,.86,.85, .84,.82,.80)
Dz22 <- c(-.57,-.10,.12,.06,-.12,.09,.05,-.01,.02,.03)
tmp <- data.frame(acf=c(1,z22, 1,Dz22),
                  lag=c(0:10, 0:10),
                  type=factor(rep(c("acf(Z)","acf(diff(Z))"), c(11,11)),
                    levels=c("acf(Z)","acf(diff(Z))")))
## hhpdf("co697ts2.pdf", width=7, height=3.5)
xyplot(acf ~ lag | type, data=tmp,
       n.used=100,
       panel=panel.acf, lwd=2,
       par.strip.text=list(cex=1.5),
       strip=function(...) strip.default(..., style=1),
       between=list(x=1, y=1),
       scales=list(alternating=FALSE, x=list(cex=1.2), y=list(cex=1.2)),
       xlab=list(cex=1.2), ylab=list("",cex=1.2),
       main=list("Time Series Question, n=100, Z.bar=60", cex=1.6),
       layout=c(2,1))
## hhdev.off()


###################################################
### code chunk number 36: tser.tex:2436-2466
###################################################
z24    <- c(.99,.94,.87,.81,.75,.65, .55, .53,.43,.40)
Dz24   <- c(.43,.28,.51,.80,.65,.44, .31, .77,.30,.20)
D4z24  <- c(.72,.67,.55,.32,.38,.23, .24, .23,.18,.13)
DD4z24 <- c(.30,.07,.32,.50,.20,.01,-.05,-.01,.02,.03)

tmp <- data.frame(acf=c(1,z24, 1,Dz24, 1,D4z24, 1,DD4z24),
                  lag=rep(0:10,4),
                  type=factor(rep(c(
                    "acf(Z)",
                    "acf(diff(Z))",
                    "acf(diff(Z,4))",
                    "acf(diff(diff(Z,4)))"),
                    c(11,11,11,11)),
                    levels=c("acf(Z)",
                    "acf(diff(Z))",
                    "acf(diff(Z,4))",
                    "acf(diff(diff(Z,4)))")))

## hhpdf("co697ts3.pdf", width=10, height=3.5)
xyplot(acf ~ lag | type, data=tmp,
       n.used=100,
       panel=panel.acf, lwd=2,
       par.strip.text=list(cex=1.2),
       strip=function(...) strip.default(..., style=1),
       between=list(x=1, y=1),
       scales=list(alternating=FALSE, x=list(cex=1.2), y=list(cex=1.2)),
       xlab=list(cex=1.2), ylab=list("",cex=1.2),
       main=list("Time Series Question, n=100, Z.bar=55", cex=1.6),
       layout=c(4,1))
## hhdev.off()


###################################################
### code chunk number 37: tser.tex:2533-2539
###################################################
data(employM16)
## hhpdf("employM16.pdf", height=8, width=7)
print(portrait=TRUE,
tsacfplots(employM16, main="employm16", xlab="", ylab="", lag.max=36, cex=1.1, lwd=1)
)
## hhdev.off()



## Mystery time series X
data(tser.mystery.X)
X <- tser.mystery.X

X.dataplot <- tsacfplots(X, lwd=1, pch.seq=16, cex=.7)
## trellis.device(file=hh("tser/figure/arima.sim.X.eps"), postscript, horizontal=TRUE); strip.background0()
X.dataplot
## dev.off()
## export.eps(hh("tser/figure/arima.sim.X.eps"))

X.loop <- if.R(r=arma.loop(X, order=c(2,0,2)),
               s=arma.loop(X, list(list(order=c(2,0,2)))))
X.diag <- rearrange.diag.arma.loop(diag.arma.loop(X.loop, x=X))
X.diagplot <- tsdiagplot(armas=X.loop, ts.diag=X.diag, lwd=1)
## trellis.device(file=hh("tser/figure/arima.sim.Xd.eps"), postscript, horizontal=TRUE); strip.background0()
X.diagplot
## dev.off()
## export.eps(hh("tser/figure/arima.sim.Xd.eps"))
X.loop
print(X.loop[["1","1"]])



## Mystery time series Y
data(tser.mystery.Y)
Y <- tser.mystery.Y

Y.dataplot <- tsacfplots(Y, lwd=1, pch.seq=16, cex=.7)
## trellis.device(file=hh("tser/figure/arima.sim.Y.eps"), postscript, horizontal=TRUE); strip.background0()
Y.dataplot
## dev.off()
## export.eps(hh("tser/figure/arima.sim.Y.eps"))

Y.loop <- if.R(r=arma.loop(Y, order=c(2,0,2)),
               s=arma.loop(Y, list(list(order=c(2,0,2)))))
Y.diag <- rearrange.diag.arma.loop(diag.arma.loop(Y.loop, x=Y))
Y.diagplot <- tsdiagplot(armas=Y.loop, ts.diag=Y.diag, lwd=1)
## trellis.device(file=hh("tser/figure/arima.sim.Yd.eps"), postscript, horizontal=TRUE); strip.background0()
Y.diagplot
## dev.off()
## export.eps(hh("tser/figure/arima.sim.Yd.eps"))
Y.loop
print(Y.loop[["1","1"]])



## Mystery time series Z
data(tser.mystery.Z)
Z <- tser.mystery.Z

Z.dataplot <- tsacfplots(Z, lwd=1, pch.seq=16, cex=.7)
## trellis.device(file=hh("tser/figure/arima.sim.Z.eps"), postscript, horizontal=TRUE); strip.background0()
Z.dataplot
## dev.off()
## export.eps(hh("tser/figure/arima.sim.Z.eps"))

Z.loop <- if.R(r=arma.loop(Z, order=c(2,0,2)),
               s=arma.loop(Z, list(list(order=c(2,0,2)))))
Z.diag <- rearrange.diag.arma.loop(diag.arma.loop(Z.loop, x=Z))
Z.diagplot <- tsdiagplot(armas=Z.loop, ts.diag=Z.diag, lwd=1)
## trellis.device(file=hh("tser/figure/arima.sim.Zd.eps"), postscript, horizontal=TRUE); strip.background0()
Z.diagplot
## dev.off()
## export.eps(hh("tser/figure/arima.sim.Zd.eps"))
Z.loop
print(Z.loop[["1","1"]])




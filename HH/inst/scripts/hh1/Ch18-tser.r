## product.s
## tsamstat.s
## arima.sim.XYZ.s
## nottem.s
## tsq.s
## ozone.s
## co697ts.s
## employM16.s
## splus.library/ARIMA-trellis.s
## -rwx------+   1 rmh None    90145 2005-05-31 22:14 tser/tser-060904.tex
## my.pacf.s



## product.s
data(product)

product.dataplot <- 
  tsacfplots(product, main="")
## trellis.device(file=hh("tser/figure/prodfig1.eps"), postscript, horizontal=TRUE); strip.background0()
print(product.dataplot)
## dev.off()
## export.eps(hh("tser/figure/prodfig1.eps"))

product.diffplot <- 
  tsacfplots(diff(product), main="")
## trellis.device(file=hh("tser/figure/prodfig2.eps"), postscript, horizontal=TRUE); strip.background0()
print(product.diffplot)
## dev.off()
## export.eps(hh("tser/figure/prodfig2.eps"))

product.loop <- if.R(r=arma.loop(product, order=c(2,1,2)),
                     s=arma.loop(product, list(list(order=c(2,1,2)))))
cat("Table 3.\n")
product.loop

product.diags <- diag.arma.loop(product.loop, x=product, lag.max=60)
product.diagplot <-   ## strip labels in R need 14in wide window.
tsdiagplot(armas=product.loop, diags=product.diags,
           x=product,
	   lag.lim=c(-2,60),
	   lag.x.at=seq(0,52,13))
## trellis.device(file=hh("tser/figure/prodfig4.eps"), postscript, horizontal=TRUE); strip.background0()
print(product.diagplot)
## dev.off()
## export.eps(hh("tser/figure/prodfig4.eps"))

cat("Table 5.\n")
print(product.loop[["2","2"]])

print(product.loop[["1","2"]])
print(product.loop[["2","1"]])
print(product.loop[["1","1"]])
## at which point there is no need for any terms


## try again with undifferenced data
product1.loop <- if.R(r=arma.loop(product, order=c(2,0,2)),
                     s=arma.loop(product, list(list(order=c(2,0,2)))))
cat("Table 6.\n")
product1.loop

product1.diags <- diag.arma.loop(product1.loop, lag.max=60, x=product)
product1.diagplot <- 
tsdiagplot(armas=product1.loop, diags=product1.diags,
	   lag.lim=c(-2,60),
	   lag.x.at=seq(0,52,13))
## trellis.device(file=hh("tser/figure/prodfig7.eps"), postscript, horizontal=TRUE); strip.background0()
print(product1.diagplot)
## dev.off()
## export.eps(hh("tser/figure/prodfig7.eps"))

cat("Table 8.\n")
print(product1.loop[["1","0"]])




## tsamstat.s
## Examples using Time Series Plotting Functions

## Displays for Direct Comparison of ARIMA Models
##
## The American Statistician, May 2002, Vol. 56, No. 2, pp. 131-138
## Richard M. Heiberger, Temple University
## Paulo Teles, Faculdade de Economia do Porto

## Many comments, marked HT, in this file refer to figure numbers in
## the Heiberger and Teles article.  Main titles for figures in the
## article have been suppressed for HH, but are left here as comments.

## For work on the computer screen, we use the residplot (called by
## tsdiagplot) and seqplot (called by tsacfplots) default lwd=0 and
## then export the graph to PostScript with the export.eps() command.
## For good quality reproduction for the book, we needed to use the
## optional lwd=1 and write directly to the PostScript file.

## co2 in S-Plus is monthly data from 1959 through 1990
## co2 in R is monthly data from 1959 through 1997
## The two co2 datasets are similar over the overlapping time periods

data(co2) ## Based on the builtin R data set

##HT New 1 1a 2
##HT	   main="Figure 1.  co2"
## trellis.device(file=hh("tser/figure/tsamsta1.ps"), postscript, horizontal=TRUE); strip.background0()
tsamsta1 <-
tsacfplots(co2,                   xlab=NULL, ylab=list("co2", cex=1.5),                   lag.max=36, cex=.5, main="", lwd=1)   ## black and white
print(tsamsta1)
## dev.off()
## export.eps(hh("tser/figure/tsamsta1.ps"))
tsamsta1color <-
tsacfplots(co2,                   xlab=NULL, ylab=list("co2", cex=1.5),                   lag.max=36, cex=.5, main=""        )  ## uses default lwd=0 on color GraphSheet
print(tsamsta1color)
## export.eps(hh("tser/figure/tsamsta1.color.ps"))

##HT main="Figure 1a (not in printed paper).  diff(co2,1)"
## trellis.device(file=hh("tser/figure/tsamsta1a.ps"), postscript, horizontal=TRUE); strip.background0()
tsamsta1a <-
tsacfplots(diff(co2,1),           xlab=NULL, ylab=list("diff(co2,1)", cex=1.5),           lag.max=36, cex=.5, main="", lwd=1)   ## black and white
print(tsamsta1a)
## dev.off()
## export.eps(hh("tser/figure/tsamsta1a.ps"))
tsamsta1acolor <-
tsacfplots(diff(co2,1),           xlab=NULL, ylab=list("diff(co2,1)", cex=1.5),           lag.max=36, cex=.5, main=""        )  ## uses default lwd=0 on color GraphSheet
print(tsamsta1acolor)
## export.eps(hh("tser/figure/tsamsta1a.color.ps"))

##HT main="Figure 2.  diff(diff(co2,1), 12)"
## trellis.device(file=hh("tser/figure/tsamsta2.ps"), postscript, horizontal=TRUE); strip.background0()
tsamsta2 <-
tsacfplots(diff(diff(co2,1), 12), xlab=NULL, ylab=list("diff(diff(co2,1), 12)", cex=1.5), lag.max=36, cex=.5, main="", lwd=1)   ## black and white
print(tsamsta2)
## dev.off()
## export.eps(hh("tser/figure/tsamsta2.ps"))
tsamsta2color <-
tsacfplots(diff(diff(co2,1), 12), xlab=NULL, ylab=list("diff(diff(co2,1), 12)", cex=1.5), lag.max=36, cex=.5, main=""        )  ## uses default lwd=0 on color GraphSheet
print(tsamsta2color)
## export.eps(hh("tser/figure/tsamsta2.color.ps"))


###HT lag.0=FALSE
##HT main="Figure 2alt (not in printed paper).  diff(diff(co2,1), 12)  with lag.0=FALSE"
tsamsta2a <-
tsacfplots(diff(diff(co2,1), 12), xlab=NULL, ylab=list("diff(diff(co2,1), 12)", cex=1.5), lag.max=36, cex=.5, main="", lwd=1, lag.0=FALSE)
print(tsamsta2a)

##HT table 0 not in paper
cat("Table 0\n")
ddco2.loop <- if.R(s=
                   arma.loop(co2,
                             list(list(order=c(2,1,2)),
                                  list(order=c(0,1,0), period=12)))
                   ,r=
                   arma.loop(co2,
                             order=c(2,1,2),
                             seasonal=list(order=c(0,1,0), period=12))
                   )
print(ddco2.loop)


##HT 3
##HT main="Figure 3"
ddco2.diags <- diag.arma.loop(ddco2.loop, co2, lag.max=36)
ddco2.diagplot <- 
tsdiagplot(armas=ddco2.loop, diags=ddco2.diags,
	   lag.lim=c(-2,38),
	   lag.x.at=seq(0,36,6),
	   lag.x.labels=c(0,"",12,"",24,"",36),
           main="", lwd=1)
## trellis.device(file=hh("tser/figure/tsamsta3.ps"), postscript, horizontal=TRUE); strip.background0()
print(ddco2.diagplot)   ## strip labels in R need 14in wide window.
## dev.off()
## export.eps(hh("tser/figure/tsamsta3.ps"))


##HT simplest specification of Figure 4.
## this has lots of overprinting of the axis labels.
if.R(s=
     tsdiagplot(co2, model=
                list(list(order=c(2,1,2)),
                     list(order=c(0,1,1), period=12)),
                main="Figure 4simple (not in printed paper)  with defaults")
     ,r=
     tsdiagplot(co2,   ## strip labels in R need 14in wide window.
                order=c(2,1,2),
                seasonal=list(order=c(0,1,1), period=12),
                main="Figure 4simple (not in printed paper)  with defaults")
     )

## table 1
ddco2.loopPQ <-
  if.R(s=
       arma.loop(co2, list(list(order=c(2,1,2)),
                           list(order=c(0,1,1), period=12)))
       ,r=
       arma.loop(co2, 
                 order=c(2,1,2),
                 seasonal=list(order=c(0,1,1), period=12))
       )
cat("Table 1\n")
print(ddco2.loopPQ)

ddco2.diagsPQ <- diag.arma.loop(ddco2.loopPQ, co2, lag.max=36)
##HT main="Figure 4"
ddco2.diagplot <-
tsdiagplot(armas=ddco2.loopPQ, diags=ddco2.diagsPQ,
	   lag.lim=c(-2,38),
	   lag.x.at=seq(0,36,6),
	   lag.x.labels=c(0,"",12,"",24,"",36),
           main="", lwd=1)
print(ddco2.diagplot)  ## strip labels in R need 14in wide window.
## export.eps(hh("tser/figure/tsamsta4.ps"))

ddco2.diagplot.0 <-
tsdiagplot(armas=ddco2.loopPQ, diags=ddco2.diagsPQ,
	   lag.lim=c(-2,38),
	   lag.x.at=seq(0,36,6),
	   lag.x.labels=c(0,"",12,"",24,"",36),
           lag.0=FALSE,
           main="Figure 4alt (not in printed paper)   with lag.0=FALSE.")
print(ddco2.diagplot.0)  ## strip labels in R need 14in wide window.


## table 1a
cat("Table 1a (not in printed paper).\n")
print(ddco2.loopPQ[["1","1"]])

## table 2
cat("Table 2\n")
co2.arima <- ddco2.loopPQ[["0","1"]]
##The final model is equivalent to
##co2.arima <-
##  arima.mle(co2, model=list(list(order=c(0,1,1)),
##                            list(order=c(0,1,1), period=12)))
print(co2.arima)

##forecast
co2.forecast <-
  if.R(s=
       arima.forecast(co2, model=co2.arima$model, n=12)
       ,r=
       predict(co2.arima, n.ahead = 12, se.fit = TRUE)
       )

##5
##HT main="Figure 5  co2 --- 1990 observed, 1991 forecast + 95% CI"
co2.plot.forecast <-
  if.R(r={
    co2.last.year <- ts(co2[457:468],
                        start=time(co2)[457],
                        frequency=frequency(co2))
    seqplotForecast(co2.last.year, co2.forecast,
                    x.at=seq(1997,1999,.5),
                    x.labels=c(1997,"",1998,"",1999),
                    xlim=c(1996.9,1999),
                    ylim=c(360, 370),
                    main="", lwd=1, cex=2)
  }, s={
    seqplotForecast(co2[373:384], co2.forecast,
                    x.at=seq(1990,1992,.5),
                    x.labels=c(1990,"",1991,"",1992),
                    xlim=c(1990,1992),
                    main="", lwd=1, cex=1.5)
  })
## trellis.device(file=hh("tser/figure/tsamsta5.ps"), postscript, horizontal=TRUE); strip.background0()
print(co2.plot.forecast)
## dev.off()
## export.eps(hh("tser/figure/tsamsta5.ps"))




## arima.sim.XYZ.s
## Mystery time series X
data(tser.mystery.X)
X <- tser.mystery.X

X.dataplot <- tsacfplots(X, lwd=1, pch.seq=16, cex=.7)
## trellis.device(file=hh("tser/figure/arima.sim.X.eps"), postscript, horizontal=TRUE); strip.background0()
print(X.dataplot)
## dev.off()
## export.eps(hh("tser/figure/arima.sim.X.eps"))

X.loop <- if.R(r=arma.loop(X, order=c(2,0,2)),
               s=arma.loop(X, list(list(order=c(2,0,2)))))
X.diag <- rearrange.diag.arma.loop(diag.arma.loop(X.loop, x=X))
X.diagplot <- tsdiagplot(armas=X.loop, ts.diag=X.diag, lwd=1)
## trellis.device(file=hh("tser/figure/arima.sim.Xd.eps"), postscript, horizontal=TRUE); strip.background0()
print(X.diagplot)
## dev.off()
## export.eps(hh("tser/figure/arima.sim.Xd.eps"))
print(X.loop)
print(X.loop[["1","1"]])



## Mystery time series Y
data(tser.mystery.Y)
Y <- tser.mystery.Y

Y.dataplot <- tsacfplots(Y, lwd=1, pch.seq=16, cex=.7)
## trellis.device(file=hh("tser/figure/arima.sim.Y.eps"), postscript, horizontal=TRUE); strip.background0()
print(Y.dataplot)
## dev.off()
## export.eps(hh("tser/figure/arima.sim.Y.eps"))

Y.loop <- if.R(r=arma.loop(Y, order=c(2,0,2)),
               s=arma.loop(Y, list(list(order=c(2,0,2)))))
Y.diag <- rearrange.diag.arma.loop(diag.arma.loop(Y.loop, x=Y))
Y.diagplot <- tsdiagplot(armas=Y.loop, ts.diag=Y.diag, lwd=1)
## trellis.device(file=hh("tser/figure/arima.sim.Yd.eps"), postscript, horizontal=TRUE); strip.background0()
print(Y.diagplot)
## dev.off()
## export.eps(hh("tser/figure/arima.sim.Yd.eps"))
print(Y.loop)
print(Y.loop[["1","1"]])



## Mystery time series Z
data(tser.mystery.Z)
Z <- tser.mystery.Z

Z.dataplot <- tsacfplots(Z, lwd=1, pch.seq=16, cex=.7)
## trellis.device(file=hh("tser/figure/arima.sim.Z.eps"), postscript, horizontal=TRUE); strip.background0()
print(Z.dataplot)
## dev.off()
## export.eps(hh("tser/figure/arima.sim.Z.eps"))

Z.loop <- if.R(r=arma.loop(Z, order=c(2,0,2)),
               s=arma.loop(Z, list(list(order=c(2,0,2)))))
Z.diag <- rearrange.diag.arma.loop(diag.arma.loop(Z.loop, x=Z))
Z.diagplot <- tsdiagplot(armas=Z.loop, ts.diag=Z.diag, lwd=1)
## trellis.device(file=hh("tser/figure/arima.sim.Zd.eps"), postscript, horizontal=TRUE); strip.background0()
print(Z.diagplot)
## dev.off()
## export.eps(hh("tser/figure/arima.sim.Zd.eps"))
print(Z.loop)
print(Z.loop[["1","1"]])




## nottem.s
## Nottingham Castle, mean monthly air temperature in degrees Fahrenheit
## January 1920 -- December 1939

library("MASS")  ## nottem is part of the MASS package

nottem.dataplot <- tsacfplots(nottem, lwd=1)
## trellis.device(file=hh("tser/figure/nottema.ps"), postscript, horizontal=TRUE); strip.background0()
print(nottem.dataplot)
## dev.off()
## export.eps(hh("tser/figure/nottema.ps"))

nottem.diff.dataplot <- tsacfplots(diff(nottem), lwd=1)  ## not displayed in book
print(nottem.diff.dataplot)

nottem.diff12.dataplot <- tsacfplots(diff(nottem, 12), lwd=1)
## trellis.device(file=hh("tser/figure/nottemb.ps"), postscript, horizontal=TRUE); strip.background0()
print(nottem.diff12.dataplot)
## dev.off()
## export.eps(hh("tser/figure/nottemb.ps"))


nottem.loop <- if.R(s=
                    arma.loop(nottem,  list(list(order=c(2,0,2)),
                                       list(order=c(2,1,0), period=12)))
                    ,r=
                    arma.loop(nottem,  order=c(2,0,2),
                                       seasonal=list(order=c(2,1,0), period=12),
                              method="ML")
                    )
print(nottem.loop, digits=4)

nottem.diag <- rearrange.diag.arma.loop(diag.arma.loop(nottem.loop, nottem))
nottem.diagplot <- tsdiagplot(armas=nottem.loop, ts.diag=nottem.diag, lwd=1)
## trellis.device(file=hh("tser/figure/nottemc.ps"), postscript, horizontal=TRUE); strip.background0()
print(nottem.diagplot)
## dev.off()
## export.eps(hh("tser/figure/nottemc.ps"))

print(nottem.loop[["1","0"]])

## if.R(s=detach("MASS"), r=detach("package:MASS"))



## tsq.s
data(tsq)

tsq.plot <- tsacfplots(tsq, main="tsq 1", lwd=1)
## trellis.device(file=hh("tser/figure/tsq1.ps"), postscript, horizontal=TRUE); strip.background0()
print(tsq.plot)
## dev.off()
## export.eps(hh("tser/figure/tsq1.ps"))

tmp <- tsq.plot$acf.plots
acfs <- if.R(r={
               acfs <- t(sapply(tmp$panel.args, `[[`, "y"))
               dimnames(acfs) <- list(c("acf","pacf"), 0:(dim(acfs)[2]-1))
               acfs
             }, s=
             t(matrix(tmp$y, ncol=2,
                      dimnames=list(matrix(tmp$x, ncol=2)[,1], c("acf","pacf"))))
             )
print(round(acfs[,1:6], digits=3))

tsq.loop <- if.R(s=
               arma.loop(tsq, model=list(order=c(2,0,2)))
               ,r=
               arma.loop(tsq, order=c(2,0,2))
               )

print(tsq.loop)

tsq.diag <-
  rearrange.diag.arma.loop(diag.arma.loop(tsq.loop, tsq))
tsq.diagplot <- tsdiagplot(armas=tsq.loop, ts.diag=tsq.diag, lwd=1)
## trellis.device(file=hh("tser/figure/tsq2.ps"), postscript, horizontal=TRUE); strip.background0()
print(tsq.diagplot)
## dev.off()
## export.eps(hh("tser/figure/tsq2.ps"))




## ozone.s
data(ozone)
arosa <- ozone

print(seqplot(arosa))   ## all data 1/26 -- 11/53

##missing data prevents acfplot
## 9/31 -- 11/53 has no missing data.

arosa.subset <- ts(arosa[68:335], start=c(1931,8), freq=12)
arosa.dataplot <-
  tsacfplots(arosa.subset, lwd=1,  ## 8/31 -- 11/53 has no missing data.
             ylab="Thickness in Dobson units",
             main="Thickness of the ozone layer at arosa")
## trellis.device(file=hh("tser/figure/ex0304-4.ps"), postscript, horizontal=TRUE); strip.background0()
print(arosa.dataplot)
## dev.off()
## export.eps(hh("tser/figure/ex0304-4.ps"))




## co697ts.s
##1
##Time Series Question (a)
##n=100
##zbar=25
##rss <- 256
##z98 <- 24
##z99 <- 26
##z100<-25
##forecast z101 and z102 and their 95% forecast intervals.


acf.1  <- c(.8,.61,.47,.40,.31,.21,.18,.11,.06,.01)
pacf.1 <- c(.8,.08,.00,-.11,.00,-.12,.07,.05,.01,.02)

tmp <- data.frame(acf=c(1,acf.1, 1,pacf.1),
                  lag=c(0:10,0:10),
                  type=rep(c("acf","pacf"),c(11,11)))


print(position=c(0,.3, 1,1),
xyplot(acf ~ lag | type, data=tmp,
       n.used=100,
       panel=panel.acf,
       par.strip.text=list(cex=1.5),
       strip=function(...) strip.default(..., style=1),
       between=list(x=1, y=1),
       scales=list(alternating=FALSE, x=list(cex=1.2), y=list(cex=1.2)),
       xlab=list(cex=1.2), ylab=list("",cex=1.2),
       main=list("Time Series Question, n=100, Z.bar=25", cex=1.6),
       layout=c(2,1))
)
## export.eps(hh("tser/figure/co697ts1.eps"))



##A2: n=100, zbar=60
##Time Series Question (b)
z22  <- c( .93, .92,.90,.90, .87,.86,.85, .84,.82,.80)
Dz22 <- c(-.57,-.10,.12,.06,-.12,.09,.05,-.01,.02,.03)

tmp <- data.frame(acf=c(1,z22, 1,Dz22),
                  lag=c(0:10,0:10),
                  type=factor(rep(c("acf(Z)","acf(diff(Z))"), c(11,11)),
                    levels=c("acf(Z)","acf(diff(Z))")))


print(position=c(0,.3, 1,1),
xyplot(acf ~ lag | type, data=tmp,
       n.used=100,
       panel=panel.acf,
       par.strip.text=list(cex=1.5),
       strip=function(...) strip.default(..., style=1),
       between=list(x=1, y=1),
       scales=list(alternating=FALSE, x=list(cex=1.2), y=list(cex=1.2)),
       xlab=list(cex=1.2), ylab=list("",cex=1.2),
       main=list("Time Series Question, n=100, Z.bar=60", cex=1.6),
       layout=c(2,1))
)
## export.eps(hh("tser/figure/co697ts2.eps"))





#A4: n=100, zbar=55
#Time Series Question (c)
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

print(position=c(0,.4, 1,1),
xyplot(acf ~ lag | type, data=tmp,
       n.used=100,
       panel=panel.acf,
       par.strip.text=list(cex=1.2),
       strip=function(...) strip.default(..., style=1),
       between=list(x=1, y=1),
       scales=list(alternating=FALSE, x=list(cex=1.2), y=list(cex=1.2)),
       xlab=list(cex=1.2), ylab=list("",cex=1.2),
       main=list("Time Series Question, n=100, Z.bar=55", cex=1.6),
       layout=c(4,1))
)
## export.eps(hh("tser/figure/co697ts3.eps"))




## employM16.s
## Employment data from DATA by Andrews and Herzberg.

## Table 65.1 United States of America Monthly Employment Figures 
##             for Males Aged 16-19 Years from 1948-1981 



data(employM16)
m16 <- employM16

m16.dataplot <-
  tsacfplots(m16, main="m16", xlab="", ylab="", lag.max=36, cex=.5, lwd=1)
## trellis.device(file=hh("tser/figure/employM16.eps"), postscript, horizontal=TRUE); strip.background0()
print(m16.dataplot)
## dev.off()
## export.eps(hh("tser/figure/employM16.eps"))




## splus.library/ARIMA-trellis.s
## These functions are included in the HH package



## my.pacf.s
## This function illustrates the definition of the pacf.
## Do NOT use in actual calculations!
##
## my.pacf requires a detrended series, otherwise the answer is
## nonsense as it starts losing precision after the first few lags.
##
## Old restriction in S-Plus in 2003.  It seems not to be relevant in 2012.
##    The argument z must have a class of "rts", "cts", or "its".
##    The function will not work with class "ts".

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

my.pacf(arosa.subset)
if.R(r=pacf(arosa.subset, plot=FALSE)$acf[2],
     s=acf(arosa.subset, type="partial", plot=FALSE)$acf[2])

### R code from vignette source 'TTAVignette.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: somecode (eval = FALSE)
###################################################
## acf(z, lag.max = ((nrow(TimeSerie))/2), na.action = na.pass)


###################################################
### code chunk number 2: somecode (eval = FALSE)
###################################################
## stl(z2, s.window="periodic", t.window=(F*10), na.action=na.fail)


###################################################
### code chunk number 3: somecode (eval = FALSE)
###################################################
##  Loess <- loess(param ~ time, Regularised.data, family="gaussian",
##                 span=0.25, control = loess.control(surface = "direct"),
##                 na.action=na.exclude)
##  tsLoess <- ts(predict(Loess),
##                start=(min(Regularised.data$YEARS)), deltat=freq)



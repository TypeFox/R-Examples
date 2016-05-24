library("FitAR")

#Figure 1. Plot of lynx time series using plot.ts
par(cex=1.5)
plot(lynx)
par(cex=1)

#Figure 2. Plot of lynx series using TimeSeriesPlot
TimeSeriesPlot(lynx, type="o", pch=16, ylab="# pelts", main="Lynx Trappings")



#Figure 3. Trellis plot for Ninemile series
graphics.off() #clear previous graphics
data("Ninemile")
print(TimeSeriesPlot(Ninemile, SubLength=200))



#Figure 4. Partial autocorrelation plot of lynx series 
graphics.off() #clear previous graphics
PacfPlot(log(lynx))



#Figure 5. Using SelectModel to select the best subset ARz or ARp and
#          comparing BIC and UBIC subset selection.
#
graphics.off() #clear previous graphics
layout(matrix(1:4,ncol=2),respect=TRUE)
ansBICp<-SelectModel(log(lynx),lag.max=15,Criterion="BIC", ARModel="ARp", Best=3)
ansUBICp<-SelectModel(log(lynx),lag.max=15, ARModel="ARp", Best=3)
ansBICz<-SelectModel(log(lynx),lag.max=15,Criterion="BIC", ARModel="ARz", Best=3)
ansUBICz<-SelectModel(log(lynx),lag.max=15, ARModel="ARz", Best=3)
par(mfg=c(1,1))
plot(ansBICp)
par(mfg=c(2,1))
plot(ansUBICp)
par(mfg=c(1,2))
plot(ansBICz)
par(mfg=c(2,2))
plot(ansUBICz)



#Figure 6. Logged spectral density function fitted to square-root of monthly
#          sunspot series using the non-subset AR and subset ARz.
#          AIC and BIC are used for the AR while BIC and UBIC are used
#          for the ARz. Takes about 115 seconds on 3.6 GHz Pentium PC.
graphics.off() #clear previous graphics
layout(matrix(1:4,ncol=2),respect=TRUE)
z<-sqrt(sunspots)
P<-200
pAIC<-SelectModel(z, lag.max=P, ARModel="AR", Best=1, Criterion="AIC")
ARAIC<-FitAR(z, pAIC)
par(mfg=c(1,1))
sdfplot(ARAIC)
title(main="AIC Order Selection")
pBIC<-SelectModel(z, lag.max=P, ARModel="AR", Best=1, Criterion="BIC")
ARBIC<-FitAR(z, pBIC)
par(mfg=c(1,2))
sdfplot(ARBIC)
title(main="BIC Order Selection")
SunspotMonthARzBIC<-SelectModel(z,lag.max=P, ARModel="ARz", Best=1, Criterion="BIC")
ARzBIC<-FitAR(z, SunspotMonthARzBIC)
par(mfg=c(2,1))
sdfplot(ARzBIC)
title(main="BIC Subset Selection")
SunspotMonthARzUBIC<-SelectModel(z,lag.max=P, ARModel="ARz", Best=1)
ARzUBIC<-FitAR(z, SunspotMonthARzUBIC)
par(mfg=c(2,2))
sdfplot(ARzUBIC)
title(main="UBIC Subset Selection")



#Figure 7. Comparing Box-Cox analyses using FitAR and MASS
library("MASS")
graphics.off() #clear previous graphics
layout(matrix(c(1,2,1,2),ncol=2))
pvec<-c(1,2,4,10,11)
out<-FitAR(lynx, ARModel="ARp", pvec)
BoxCox(out)
PMAX<-max(pvec)
Xy <- embed(lynx, PMAX + 1)
y <- Xy[, 1]
X <- (Xy[, -1])[, pvec] #pvec != 1
outlm<-lm(y~X)
boxcox(outlm,lambda=seq(0.0,0.6,0.05))



#Figure 8
graphics.off() #clear previous graphics
BoxCox(AirPassengers) #takes about 30 sec



#Figure 9
graphics.off() #clear previous graphics
data("rivers")
BoxCox(rivers)
title(sub="Length of 141 North American Rivers")



#Figure 10
graphics.off() #clear previous graphics
data("USTobacco")
TimeSeriesPlot(USTobacco, aspect=1)



#Figure 11
graphics.off() #clear previous graphics
outUST<-arima(USTobacco, c(0,1,1))
BoxCox(outUST)



#Figure 12. Basic diagnostic plots for ARp fitted to the log lynx series
graphics.off() #clear previous graphics
out<-FitAR(log(lynx), ARModel="ARp", c(1,2,4,10,11))
plot(out, terse=TRUE)



#Figure 13. RSF plot for ARp fitted to log lynx series
graphics.off() #clear previous graphics
out<-FitAR(log(lynx), ARModel="ARp", c(1,2,4,10,11))
rfs(out)




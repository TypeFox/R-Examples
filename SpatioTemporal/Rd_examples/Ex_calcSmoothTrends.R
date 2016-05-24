##let's load some data
data(mesa.model)

##let's compute two smooth trend functions
trend <- calcSmoothTrends(mesa.model, n.basis=2)

##or with some other parameters for the splines
trend.alt <- calcSmoothTrends(mesa.model, n.basis=2, df=100)

##and study the trends
par(mfrow=c(2,1), mar=c(2.5,2.5,.5,.5))
plot(trend$trend$date, trend$trend$V1, type="l", ylab="", xlab="",
     ylim=range(c(trend$trend$V1, trend$trend$V2)))
lines(trend$trend$date, trend$trend$V2, col=2)
plot(trend.alt$trend$date, trend.alt$trend$V1,
     type="l", ylab="", xlab="",
     ylim=range(c(trend.alt$trend$V1, trend.alt$trend$V2)))
lines(trend.alt$trend$date, trend.alt$trend$V2, col=2)

##Let's exclude locations with fewer than 100 observations
IND <- names(which(table(mesa.model$obs$ID) >= 100))
##now we also compute the CV trends.
trend2 <- calcSmoothTrends(mesa.model, n.basis=2, subset=IND, cv=TRUE)

##Let's compare to the previous result
lines(trend2$trend$date, trend2$trend$V1, lty=2)
lines(trend2$trend$date, trend2$trend$V2, lty=2, col=2)

##we can also study the cross validated results to examine the
##possible variation in the estimated trends.
plot(trend$trend$date, trend2$trend$V1, type="n", ylab="", xlab="",
     ylim=range(c(trend2$trend$V1, trend2$trend$V2)))
for(i in 1:length(trend2$trend.cv)){
  lines(trend2$trend.cv[[i]]$date, trend2$trend.cv[[i]]$V1, col=1)
  lines(trend2$trend.cv[[i]]$date, trend2$trend.cv[[i]]$V2, col=2)
}

##trend functions for every week (mesa.model$obs$date is every 2-weekss)
trend.more <- calcSmoothTrends(mesa.model, n.basis=2,
                               extra.dates=seq(min(mesa.model$obs$date),
                                 max(mesa.model$obs$date), by=7))
##This results in a message detailing how many times that
##have hade interpolated trends (i.e. no data)

##compare to the earlier
plot(trend$trend$date, trend$trend$V1, pch=19,
     ylim=range(c(trend$trend$V1, trend.more$trend$V1)))
points(trend.more$trend$date, trend.more$trend$V1, col=2, pch=3)

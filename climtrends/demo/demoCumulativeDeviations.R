# example of cumulative deviations test
# Dury, G. (1980). "Step-functional changes in precipitation at Sydney", Australian Geographical Studies 18, 62-78.
# Rainfall data from Madison (Wisconsin)

data(annual.precipitation.totals.Madison)
Madison <- annual.precipitation.totals.Madison
Madison[,2]<-Madison[,2]*25.4
Madison <- Madison[-(1:6),]
Madison <- Madison[1:57,]
par(mfrow=c(2,1),mar=c(1,1,1,1),oma=c(1, 3, 0,0))
plot(Madison,type='s')
y <- cumulativeDeviations(Madison[,2])
m <- mean(Madison[,2])
n <- dim(Madison)[1]
plot(Madison[,1],y,type='p',col=2 )

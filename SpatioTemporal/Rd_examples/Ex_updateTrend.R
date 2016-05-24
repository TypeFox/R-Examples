##load data
data(mesa.model)

##default data and time trend for one location
par(mfrow=c(3,1), mar=c(2.5,2.5,3,1))
plot(mesa.model)

##let's try with no trend
mesa.model.0 <- updateTrend(mesa.model, n.basis=0)
plot(mesa.model.0)

##...and two basis functions, based on only AQS sites and much less smooth
subset <- mesa.model$locations$ID[mesa.model$locations$type=="AQS"]
mesa.model.2 <- updateTrend(mesa.model, n.basis=2, subset=subset, df=100)
plot(mesa.model.2)

##Compute trends based on only 10 sites (and compute the cross-validated
##trends leaving each of the sites out
smooth.trend <- calcSmoothTrends(mesa.model, n.basis=2, cv=TRUE,
                                 subset=subset[1:10])

##update trends using the function definition
mesa.model <- updateTrend(mesa.model, fnc=smooth.trend$trend.fnc)

##and create objects with each of the trends.
mesa.model.cv <- vector("list", length(smooth.trend$trend.fnc.cv))
for(i in 1:length(mesa.model.cv)){
  suppressMessages(mesa.model.cv[[i]] <- updateTrend(mesa.model, 
                   fnc=smooth.trend$trend.fnc.cv[[i]]))
}

##plot
par(mfrow=c(1,1),mar=c(2.5,2.5,3,1))
plot(mesa.model)
for(i in 1:length(mesa.model.cv)){
  plot(mesa.model.cv[[i]], add=TRUE, col=i, pch=NA, lty=c(NA,2))
}

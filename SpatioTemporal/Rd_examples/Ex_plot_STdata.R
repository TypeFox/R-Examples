##load data
data(mesa.model)

##default plot
plot(mesa.model)

##plot monitor locations
plot(mesa.model, "loc")

##different names/colours/etc
plot(mesa.model, "loc", main="A nice plot", col=c("green","blue"),
    legend.names=c("Sites of one type", "..and of the other"),
    legend.loc="topleft")

##composite time-trend
plot(mesa.model, "loc.obs", legend.loc="bottomleft", cex=.5, pch=c(3,4))

##plot tim-series for the first site,
layout(matrix(c(1,2,3,1,2,4),3,2))
plot(mesa.model, "obs", ID=1, col=c("red", "black"))
##residuals from the temporal trends,
plot(mesa.model, "res", ID=1, col=c("black","grey"))
##afc 
plot(mesa.model, "acf", ID=1)
##... and pafc for the residuals
plot(mesa.model, "pacf", ID=1, ci.col="red")

##Different site and with no temporal trend.
mesa.model <- updateTrend(mesa.model, n.basis=0)
layout(matrix(c(1,2,3,1,2,4),3,2))
plot(mesa.model, "obs", ID="60370016")
plot(mesa.model, "res", ID="60370016")
plot(mesa.model, "acf", ID="60370016")
plot(mesa.model, "pacf", ID="60370016")

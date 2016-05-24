##load the data
data(mesa.data.raw)
##and create STdata-object
mesa.data <- createSTdata(mesa.data.raw$obs, mesa.data.raw$X, n.basis=2,
                          SpatioTemporal=mesa.data.raw["lax.conc.1500"])

##plot time-series for the first site,
par(mfrow=c(3,2),mar=c(2.5,2.5,3,1))
plot(mesa.data, "obs", ID=1)
##And combined for all sites
plot(mesa.data, "loc.obs", legend.loc="bottomleft")

##attempt to detrend
mesa.data.detrend <- detrendSTdata(mesa.data)
##examine object, note the trends
mesa.data.detrend

##plot detrended time-series for the first site,
plot(mesa.data.detrend, "obs", ID=1)
##And combined for all sites
plot(mesa.data.detrend, "loc.obs", legend.loc="bottomleft")

##use different detrending for different types of locations
mesa.data.detrend2 <- detrendSTdata(mesa.data, region=mesa.data$covars$type)
##examine object, note the trends
mesa.data.detrend2
##plot for the first site,
plot(mesa.data.detrend2, "obs", ID=1)
plot(mesa.data.detrend2, "loc.obs", legend.loc="bottomleft")

##compare the two fitted and removed trends
print(mesa.data.detrend$fit.trend)
print(mesa.data.detrend2$fit.trend)

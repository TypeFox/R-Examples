##load the data
data(mesa.data.raw)
##and create STdata-object
mesa.data <- createSTdata(mesa.data.raw$obs, mesa.data.raw$X, n.basis=2,
                          SpatioTemporal=mesa.data.raw["lax.conc.1500"])

##specify locations, using x/y and specifying long/lat and picking
##type as an additional field
loc.spec <- list(coords=c("x","y"), long.lat=c("long","lat"), others="type")
##create the location data.frame
str( processLocation(mesa.data, loc.spec) )

##specify only locations
str( processLocation(mesa.data, list(coords=c("x","y"))) )

##different coordinates for beta and nu fields
loc.spec <- list(coords=c("x","y"), coords.nu=c("long","lat"))
str( processLocation(mesa.data, loc.spec) )

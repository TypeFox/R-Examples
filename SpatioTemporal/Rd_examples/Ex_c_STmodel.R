##load the data
data(mesa.data.raw)
##and create STdata-object
mesa.data <- createSTdata(mesa.data.raw$obs, mesa.data.raw$X, n.basis=2,
                          SpatioTemporal=mesa.data.raw["lax.conc.1500"])

##keep only observations from the AQS sites
ID.AQS <- mesa.data$covars$ID[ mesa.data$covars$type=="AQS" ]
mesa.data$obs <- mesa.data$obs[mesa.data$obs$ID %in% ID.AQS,]

##model specification
LUR <- list(~log10.m.to.a1 + s2000.pop.div.10000 + km.to.coast,
            ~km.to.coast, ~km.to.coast)
locations <- list(coords=c("x","y"), long.lat=c("long","lat"), others="type")

##create reduced model, without and with a spatio-temporal covariate.
mesa.model <- createSTmodel(mesa.data, LUR=LUR, locations=locations,
                            strip=TRUE)
mesa.model.ST <- createSTmodel(mesa.data, LUR=LUR, ST=1,
                               locations=locations, strip=TRUE)
##and non stripped version
mesa.model.full <- createSTmodel(mesa.data, LUR=LUR, ST=1,
                                 locations=locations)

##combine, this adds the missing locations
mesa.model$locations$ID
c(mesa.model, mesa.data)$locations$ID

##or we could study the summary output
print(c(mesa.model.ST, mesa.data))

##no change since we're tryin to adding existing sites
mesa.model.full$locations$ID
c(mesa.model.full, mesa.data)$locations$ID

##We can also combine two STmodels
print(c(mesa.model, mesa.model.full))

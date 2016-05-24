##load the data
data(mesa.data.raw)
##and create STdata-object
mesa.data <- createSTdata(mesa.data.raw$obs, mesa.data.raw$X, n.basis=2,
                          SpatioTemporal=mesa.data.raw["lax.conc.1500"])

##define land-use covariates
LUR <-  list(~log10.m.to.a1+s2000.pop.div.10000+km.to.coast,
             ~km.to.coast, ~km.to.coast)
##and covariance model
cov.beta <- list(covf="exp", nugget=FALSE)
cov.nu <- list(covf="exp", nugget=TRUE, random.effect=FALSE)
##which locations to use
locations <- list(coords=c("x","y"), long.lat=c("long","lat"), others="type")

##create object
mesa.model <- createSTmodel(mesa.data, LUR=LUR, ST="lax.conc.1500",
                            cov.beta=cov.beta, cov.nu=cov.nu,
                            locations=locations)
print(mesa.model)
##This is the same as data(mesa.model)

##lets try some alternatives:
model.none <- createSTmodel(mesa.data, LUR=NULL, ST=NULL)
print(model.none)

##Specify LUR:s using numbers
names(mesa.data$covars)
model.diff <- createSTmodel(mesa.data, LUR=list(c(7,10,11,12),11:12,11:12),
                            ST=1)
print(model.diff)

##Same covariates for all temporal trends, calling by name
##but with different covariance models for each trend, and nugget that depends
##on monitor type
model.same <- createSTmodel(mesa.data, LUR=c("log10.m.to.a1", "log10.m.to.road",
                                         "km.to.coast","s2000.pop.div.10000"),
                            ST="lax.conc.1500", cov.nu=list(nugget="type"),
                            cov.beta=list(covf=c("exp","exp2","iid"),
                              nugget=c(FALSE, FALSE, TRUE)) )
print(model.same)

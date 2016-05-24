##load the data
data(mesa.data.raw)
##and create STdata-object
mesa.data <- createSTdata(mesa.data.raw$obs, mesa.data.raw$X, n.basis=2,
                          SpatioTemporal=mesa.data.raw["lax.conc.1500"])

##create a simple set of covariates
processLUR(mesa.data, list(c(7:9),7,8))

##or a structure with the same covariates for all
##temporal trends
processLUR(mesa.data, c(7,11))

##or a structure with only intercept for the temporal trends
processLUR(mesa.data, list(c(7:9),NULL,NULL))

##Ask for covariates by name
processLUR(mesa.data, list(c("log10.m.to.a1","log10.m.to.a2"),
                           "log10.m.to.a1","log10.m.to.a1"))
##use formula for part of it
processLUR(mesa.data, list(~log10.m.to.a1+log10.m.to.a2+log10.m.to.a1*km.to.coast,
                           "log10.m.to.a1", "log10.m.to.a1"))

##Ask for non-existent covariate by name or formula, or location
##for each temporal trend)
try(processLUR(mesa.data, list("log10.m.to.a4",~log10.m.to.a1+log10.m.to.a4, 25)))

##create a simple set of spatio-temporal covariates
processST(mesa.data, 1)
##or create a empty set of spatio-temporal covariates
processST(mesa.data, NULL)
##by name
processST(mesa.data, "lax.conc.1500")


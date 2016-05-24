##load the raw data
data(mesa.data.raw)

##extract observations and covariates
obs <- mesa.data.raw$obs
covars <- mesa.data.raw$X

##list with the spatio-temporal covariates
ST.list <- list(lax.conc.1500=mesa.data.raw$lax.conc.1500)

##create STdata object
mesa.data <- createSTdata(obs, covars, SpatioTemporal=ST.list)
print(mesa.data)

##create object with mean 0 spatio temporal covariate
mesa.data.2 <- createSTdata(obs, covars, SpatioTemporal=ST.list,
                            mean.0.ST=TRUE)
print(mesa.data.2)

##create object with mean 0 spatio temporal covariate, and
##trend with two components, and additional dates (every seventh day)
extra.dates <- seq(min(as.Date(rownames(obs))),
                   max(as.Date(rownames(obs))), by=7)
mesa.data.3 <- createSTdata(obs, covars, n.basis=2, extra.dates=extra.dates)
print(mesa.data.3)

##load the data
data(mesa.data.raw)
##and create STdata-object
mesa.data <- createSTdata(mesa.data.raw$obs, mesa.data.raw$X, n.basis=2,
                          SpatioTemporal=mesa.data.raw["lax.conc.1500"])

mesa.data.mean0 <- removeSTcovarMean(mesa.data)

##compare the data structures
##geographic covariates
summary(mesa.data$covars)
summary(mesa.data.mean0$covars)

##mean of the spatio-temporal covariate, note that the new
##contains both mean-zero and original
cbind(colMeans(mesa.data$SpatioTemporal),
      colMeans(mesa.data.mean0$SpatioTemporal))

\dontshow{
  if( max(abs(colMeans(mesa.data.mean0$SpatioTemporal[,,2]))) > 1e-10 ){
    stop("remove.ST.mean 1: mean(centred ST) != 0")
  }
  if( max(abs(mesa.data.mean0$covars$mean.lax.conc.1500 - 
              colMeans(mesa.data$SpatioTemporal))) > 1e-10 ){
    stop("remove.ST.mean 2: extracted mean != mean(ST)")
  }
}

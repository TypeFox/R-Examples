#######  examples from User's Guide 

  require("dse")
 
 data("egJofF.1dec93.data", package="dse") 
 
  eg4.DSE.data<- egJofF.1dec93.data
  outputData(eg4.DSE.data) <- outputData(eg4.DSE.data, series=c(1,2,6,7))

  eg4.DSE.model <- estVARXls(eg4.DSE.data)
 
  new.data <- TSdata(
              input= ts(rbind(inputData(eg4.DSE.data), matrix(.1,10,1)), 
                       start=start(eg4.DSE.data),
                       frequency=frequency(eg4.DSE.data)),    
              output=ts(rbind(outputData(eg4.DSE.data),matrix(.3,5,4)), 
                       start=start(eg4.DSE.data),
                       frequency=frequency(eg4.DSE.data)))
  seriesNames(new.data) <- seriesNames(eg4.DSE.data)

  z  <- l(TSmodel(eg4.DSE.model), trimNA(new.data)) 

  cat("Forecasting...\n")

  zz <- forecast(TSmodel(eg4.DSE.model), new.data)
  z <-  forecast(TSmodel(eg4.DSE.model), trimNA(new.data), 
		conditioning.inputs=inputData(new.data))
  tfplot(zz, start=c(1990,6))

  z <- forecast(eg4.DSE.model, conditioning.inputs.forecasts=matrix(.5,6,1)) 
  summary(z)
  tfplot(z)

  tfplot(z, start=c(1990,1)) 

  z <- l(TSmodel(eg4.DSE.model), new.data)
  tfplot(z)
  
  cat("Forecasting from different points in time...\n")

  z <- featherForecasts(TSmodel(eg4.DSE.model), new.data)
  tfplot(z)
  
  zz <-featherForecasts(TSmodel(eg4.DSE.model), new.data,
                          from.periods =c(20,50,60,70,80), horizon=150)
  tfplot(zz)

  cat("Forecasting for different horizons...\n")

  z <- horizonForecasts(TSmodel(eg4.DSE.model), new.data, horizons=c(1,3,6))
  tfplot(z)

  cat("Forecasting error...\n")

  fc1 <- forecastCov(TSmodel(eg4.DSE.model), data=eg4.DSE.data)
  tfplot(fc1, mar=c(2.1,4.1,4.1,2.1))

  tfplot(forecastCov(TSmodel(eg4.DSE.model), data=eg4.DSE.data, horizons= 1:4),
     mar=c(2.1,4.1,4.1,2.1) ) 
 
  fc2 <- forecastCov(TSmodel(eg4.DSE.model), data=eg4.DSE.data, zero=T, trend=T)
  tfplot(fc2, mar=c(2.1,4.1,4.1,2.1) )


data(presidentialForecast)

full.forecasts<-presidentialForecast[,c(1:6)]
full.observed<-presidentialForecast[,7]

# Prediction for the 2008 election
this.ForecastData<-makeForecastData(.predCalibration=full.forecasts[c(1:14),],.outcomeCalibration=full.observed[c(1:14)],.predTest=full.forecasts[15,], .outcomeTest=full.observed[15], .modelNames=c("Campbell", "Lewis-Beck","EWT2C2","Fair","Hibbs","Abramowitz"))

thisEnsemble<-calibrateEnsemble(this.ForecastData, model="normal", useModelParams=FALSE, tol=0.000000001)

summary(thisEnsemble, period="calibration", showCoefs=FALSE)
plot(thisEnsemble, period="test", main="Forecast of 2008 Election")



 

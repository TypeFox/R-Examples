data(calibrationSample)
data(testSample)

this.ForecastData <- makeForecastData(.predCalibration=calibrationSample[,c("LMER", "SAE", "GLM")],
                                      .outcomeCalibration=calibrationSample[,"Insurgency"],
                                      .predTest=testSample[,c("LMER", "SAE", "GLM")],
                                      .outcomeTest=testSample[,"Insurgency"],
                                      .modelNames=c("LMER", "SAE", "GLM"))

this.ensemble <- calibrateEnsemble(this.ForecastData, model="logit", tol=0.0001, maxIter=25000, exp=3)

summary(this.ensemble, period="calibration")
plot(this.ensemble, period="calibration")
summary(this.ensemble, period="test", showCoefs=FALSE)
plot(this.ensemble, period="test")




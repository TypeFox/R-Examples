## ----echo = FALSE--------------------------------------------------------
library(weathermetrics)

## ------------------------------------------------------------------------
data(lyon)
lyon$T_F <- celsius.to.fahrenheit(lyon$TemperatureC)
lyon$DP_F <- celsius.to.fahrenheit(lyon$DewpointC)
lyon

## ------------------------------------------------------------------------
data(norfolk)
norfolk$T_C <- fahrenheit.to.celsius(norfolk$TemperatureF)
norfolk$DP_C <- fahrenheit.to.celsius(norfolk$DewpointF)
norfolk

## ------------------------------------------------------------------------
data(lyon)
lyon$RH <- dewpoint.to.humidity(t = lyon$TemperatureC,
                                dp = lyon$DewpointC,
                                temperature.metric = "celsius")
lyon

## ------------------------------------------------------------------------
data(newhaven)
newhaven$DP <- humidity.to.dewpoint(t = newhaven$TemperatureF,
                                    rh = newhaven$Relative.Humidity,
                                    temperature.metric = "fahrenheit")
newhaven

## ------------------------------------------------------------------------
data(newhaven)
newhaven$DP <- humidity.to.dewpoint(t = newhaven$TemperatureF,
                                    rh = newhaven$Relative.Humidity,
                                    temperature.metric = "fahrenheit")
newhaven$DP_C <- fahrenheit.to.celsius(newhaven$DP)
newhaven

## ------------------------------------------------------------------------
data(suffolk)
suffolk$HI <- heat.index(t = suffolk$TemperatureF,
                         rh = suffolk$Relative.Humidity,
                         temperature.metric = "fahrenheit",
                         output.metric = "fahrenheit")
suffolk

## ------------------------------------------------------------------------
data(lyon)
lyon$HI_F <- heat.index(t = lyon$TemperatureC,
                      dp = lyon$DewpointC,
                      temperature.metric = "celsius",
                      output.metric = "fahrenheit")
lyon

## ------------------------------------------------------------------------
df <- data.frame(T = c(NA, 90, 85),
                 DP = c(80, NA, 70))
df$RH <- dewpoint.to.humidity(t = df$T, dp = df$DP,
                              temperature.metric = "fahrenheit")
df

## ------------------------------------------------------------------------
df <- data.frame(T = c(90, 90, 85),
                 DP = c(80, 95, 70))
df$heat.index <- heat.index(t = df$T, dp = df$DP,
                            temperature.metric = 'fahrenheit')
df

## ------------------------------------------------------------------------
data(suffolk)
suffolk$TempC <- fahrenheit.to.celsius(suffolk$TemperatureF,
                                       round = 5)
suffolk$HI <- heat.index(t = suffolk$TemperatureF, 
                         rh = suffolk$Relative.Humidity,
                         round = 3)
suffolk


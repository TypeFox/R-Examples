## ---- message=FALSE------------------------------------------------------
library(water)


## ------------------------------------------------------------------------
aoi <- createAoi(topleft = c(273110, -3914450), 
                 bottomright = c( 288050, -3926650), EPSG = 32619)

## ------------------------------------------------------------------------
csvfile <- system.file("extdata", "apples.csv", package="water")
MTLfile <- system.file("extdata", "L7.MTL.txt", package="water")
WeatherStation <- read.WSdata(WSdata = csvfile, date.format = "%d/%m/%Y",
                  lat=-35.42222, long= -71.38639, elev=201, height= 2.2,
                  MTL = MTLfile)

## ---- fig.width = 5------------------------------------------------------
print(WeatherStation)

plot(WeatherStation, alldata=FALSE)

## ---- fig.width = 5------------------------------------------------------
image.DN <- L7_Talca[[c(1:5,7)]]
B6 <- L7_Talca[[6]]
plot(image.DN)

## ---- fig.width = 5, warning=FALSE---------------------------------------
Energy.Balance <- METRIC.EB(image.DN = image.DN, plain=TRUE, 
                            WeatherStation = WeatherStation, 
                            ETp.coef = 1.2, MTL=MTLfile, 
                            sat="L7", thermalband=B6)

## ---- fig.width = 5------------------------------------------------------
plot(Energy.Balance)

## ------------------------------------------------------------------------
ET_WS <- dailyET(WeatherStation = WeatherStation)

## ---- fig.width = 5, warning=FALSE---------------------------------------
ET.24 <- ET24h(Rn=Energy.Balance$NetRadiation, G=Energy.Balance$SoilHeat, 
               H=Energy.Balance$SensibleHeat, 
               Ts=Energy.Balance$surfaceTemperature, 
               WeatherStation = WeatherStation, ETr.daily=ET_WS)


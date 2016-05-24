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

## ------------------------------------------------------------------------
checkSRTMgrids(image.DN)

## ------------------------------------------------------------------------
DEM <- DEM_Talca

## ---- fig.width = 5------------------------------------------------------
surface.model <-METRICtopo(DEM)

solar.angles.r <- solarAngles(surface.model = surface.model, 
                              WeatherStation = WeatherStation, MTL = MTLfile)

plot(solar.angles.r)

Rs.inc <- incSWradiation(surface.model = surface.model, 
                         solar.angles = solar.angles.r, 
                         WeatherStation = WeatherStation)

## ---- fig.width = 5------------------------------------------------------
image.TOAr <- calcTOAr(image.DN = image.DN, sat="L7", MTL = MTLfile, 
                       incidence.rel = solar.angles.r$incidence.rel)

image.SR <- calcSR(image.TOAr=image.TOAr, sat = "L7", 
                   surface.model=surface.model, 
                   incidence.hor = solar.angles.r$incidence.hor, 
                   WeatherStation=WeatherStation, ESPA = F)

albedo <- albedo(image.SR = image.SR,  coeff="Tasumi", sat = "L7")

## ---- fig.width = 5------------------------------------------------------
LAI <- LAI(method = "metric2010", image = image.TOAr, L=0.1, sat = "L7")

plot(LAI)

## ---- warning=FALSE, fig.width = 5---------------------------------------
Ts <- surfaceTemperature(LAI=LAI, sat = "L7", thermalband = B6,
                         WeatherStation = WeatherStation)

Rl.out <- outLWradiation(LAI = LAI, Ts=Ts)

Rl.inc <- incLWradiation(WeatherStation,DEM = surface.model$DEM, 
                         solar.angles = solar.angles.r, Ts= Ts)

## ---- fig.width = 5------------------------------------------------------
Rn <- netRadiation(LAI, albedo, Rs.inc, Rl.inc, Rl.out)

plot(Rn)

## ---- fig.width = 5------------------------------------------------------
G <- soilHeatFlux(image = image.SR, Ts=Ts,albedo=albedo, 
                  Rn=Rn, LAI=LAI, sat = "L7")

plot(G)

## ---- fig.width = 5------------------------------------------------------
Z.om <- momentumRoughnessLength(LAI=LAI, mountainous = TRUE, 
                                method = "short.crops", 
                                surface.model = surface.model)

hot.and.cold <- calcAnchors(image = image.TOAr, Ts = Ts, LAI = LAI, plots = F,
                            albedo = albedo, Z.om = Z.om, n = 1, 
                            anchors.method = "CITRA-MCB", sat = "L7", 
                            deltaTemp = 5, verbose = FALSE)

H <- calcH(anchors = hot.and.cold, Ts = Ts, Z.om = Z.om, 
           WeatherStation = WeatherStation, ETp.coef = 1.05, sat = "L7", 
           Z.om.ws = 0.0018, DEM = DEM, Rn = Rn, G = G, verbose = FALSE)

## ------------------------------------------------------------------------
ET_WS <- dailyET(WeatherStation = WeatherStation)

## ---- fig.width = 5------------------------------------------------------
ET.24 <- ET24h(Rn, G, H$H, Ts, WeatherStation = WeatherStation, ETr.daily=ET_WS)


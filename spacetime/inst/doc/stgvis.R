## ----setOptions, message=FALSE-------------------------------------------
library(googleVis)
op <- options(gvis.plot.tag='chart')

## ----results='asis', eval=TRUE-------------------------------------------
library(spacetime)
data(air) # loads rural and DE_NUTS1
rural = STFDF(stations, dates, data.frame(PM10 = as.vector(air)))
library(ISOcodes)
data("ISO_3166_2")
## State names are already in German
ISO_3166_2_DE <- subset(ISO_3166_2, Country %in% "DE")
plot(
  gvisTable(ISO_3166_2_DE)
  )

## ----results='asis'------------------------------------------------------
ISO_3166_2_DE <- ISO_3166_2_DE[order(ISO_3166_2_DE$Name),]
DE_NUTS1$name <- ISO_3166_2_DE$Name

## ----GeoMapExample, results='asis', tidy=FALSE---------------------------
library(googleVis)
## Create list with options for Geo Chart to be used
geoChartDE <- list(region="DE", 
                   resolution="provinces",
                   legend="{numberFormat:'#,###.00'}") 
plot(
  gvisGeoChart(DE_NUTS1@data, locationvar = "name", 
               colorvar = "Shape_Area",
               options=geoChartDE)
  )

## ----results='asis'------------------------------------------------------
DE_NUTS1$code = ISO_3166_2_DE$Code
plot(
  gvisGeoChart(DE_NUTS1@data, locationvar = "code", 
               colorvar = "Shape_Area",
               options=geoChartDE)
  )

## ----results='asis'------------------------------------------------------
DE_NUTS1.years = STF(DE_NUTS1, as.Date(c("2008-01-01", "2009-01-01")))
agg = aggregate(rural[,"2008::2009"], DE_NUTS1.years, mean, na.rm=TRUE)
d = agg[,1]@data # select time step one, take attr table of resulting SpatialPolygonsDataFrame object
d$code = ISO_3166_2_DE$Code # add region codes
d$PM10 = round(d$PM10, 2) # avoid 12-digit numbers printed when hovering over

M <- gvisGeoChart(na.omit(d), locationvar = "code", 
                 colorvar = "PM10",
                 options=c(geoChartDE, height=400)) # drop NA values for Geo chart

## Add German state names for table
d$code <- ISO_3166_2_DE$Name
Tbl <- gvisTable(d, options=list(height=400), 
                 formats=list(PM10="#,###.0"))
plot(gvisMerge(M, Tbl, horizontal=TRUE))

## ----results='asis'------------------------------------------------------
library(sp)
data("wind", package = "gstat")
wind.loc$y = as.numeric(char2dms(as.character(wind.loc[["Latitude"]])))
wind.loc$x = as.numeric(char2dms(as.character(wind.loc[["Longitude"]])))
wind.loc$yx = paste(wind.loc$y, wind.loc$x, sep=":")
m = round(apply(wind,2,mean), 3)
sd = round(apply(wind,2,sd), 3)
wind.loc$mean = m[match(wind.loc$Code, names(m))]
wind.loc$sd = sd[match(wind.loc$Code, names(sd))]
plot( gvisGeoChart(wind.loc, "yx", "mean", "sd", "Station",
                   options = list(region = "IE", 
                                  legend="{numberFormat:'#.00'}"))
      )

## ----results='asis'------------------------------------------------------
# select:
wind = wind[wind$year > 76,] # select three years, 1976-1978
time = ISOdate(wind$year+1900, wind$month, wind$day)
wind.stack = stack(wind[,4:15])
names(wind.stack) = c("wind_speed", "station")
wind.stack$time = rep(time, 12)
plot(
  gvisAnnotatedTimeLine(wind.stack, "time", "wind_speed", "station",
                        options = list(width = "1000px", 
                                       height = "500px",
                                       zoomStartTime=time[length(time)-365], 
                                       zoomEndTime=max(time)))
)

## ----results='asis'------------------------------------------------------
wind$time = time
plot(
  gvisLineChart(wind[1:200,], "time", names(wind)[4:9])
  )

## ----results='asis'------------------------------------------------------
wind$time = time
wind[4:5, 7:9] = NA
plot(
  gvisLineChart(wind[1:10,], "time", names(wind)[4:9])
  )

## ----results='asis'------------------------------------------------------
library(spacetime)
data(air)
d = as(rural[1:10,"2001"], "data.frame") # daily, 2001
plot(
  gvisAnnotationChart(d, "time", "PM10", "sp.ID",
                      options = list(width = "1000px", 
                                     height = "500px"))
  # zoomStartTime=d$time[1], zoomEndTime=d$time[1]+365))
  )

## ----MotionChartExample, results='asis', tidy=FALSE----------------------
data("Produc", package = "plm")
plot(
  gvisMotionChart(Produc, idvar="state", timevar="year")
)

## ----resetOptions--------------------------------------------------------
## Set options back to original options
options(op)


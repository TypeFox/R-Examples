library(RgoogleMaps)
data(eq2010, package = 'MSG')
boundingbox = qbbox(lon = eq2010[ , "long"], lat = eq2010[ , "lat"])
EQMap = GetMap.bbox(boundingbox$lonR, boundingbox$latR, destfile = "EQ.png", maptype = "satellite", zoom = 6)
PlotOnStaticMap(EQMap, lon = eq2010$long, lat = eq2010$lat, pch = 20, col = 'red',
                verbose = 0)

bgcol = rgb(255, 0, 0, alpha = 95, maxColorValue = 255)
PlotOnStaticMap(EQMap, lon = eq2010$lon, lat = eq2010$lat, pch = 21, bg = bgcol,
                col = "darkred", cex = 1.25 * eq2010$ms, verbose = 0)

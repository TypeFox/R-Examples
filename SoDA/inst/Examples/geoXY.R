"geoXY" <-
function(latitude, longitude,
         lat0 = min(latitude, na.rm=TRUE),
         lon0 = min(longitude, na.rm=TRUE),
         unit = 1.) {
    lat0 <- rep(lat0, length(latitude))
    lon0 <- rep(lon0, length(longitude))
    yDist <- geoDist(lat0, lon0, latitude, lon0)
    xDist <- geoDist(lat0, lon0, lat0, longitude)
    cbind(X= xDist * sign(longitude - lon0)/unit,
              Y = yDist * sign(latitude - lat0)/unit)
}

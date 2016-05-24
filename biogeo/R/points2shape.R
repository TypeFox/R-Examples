points2shape <-
function (dat, x, y, fn) 
{
    cn <- names(dat)
    f1 <- match(x, cn)
    f2 <- match(y, cn)
    xn <- coord2numeric(dat[, f1])
    yn <- coord2numeric(dat[, f2])
    coords <- cbind(xn, yn)
    data <- as.data.frame(dat)
    z1 <- SpatialPointsDataFrame(coords, data, bbox = bbox(coords), 
        proj4string = CRS("+proj=longlat +datum=WGS84"))
    writePointsShape(z1, fn)
}

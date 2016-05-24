`ensemble.area` <- function(
    x=NULL, km2=TRUE
)
{
    #    if (! require(dismo)) {stop("Please install the dismo package")}
    if(inherits(x, "RasterLayer") == F) {stop("x is not a RasterLayer object")}
    if(raster::isLonLat(x) == F) {stop("x is not in longitude-latitude coordinates")}
    cat(paste("\n", "Cell frequencies", "\n", sep = ""))
    print(raster::freq(x))
    count.polygon <- raster::rasterToPolygons(x, dissolve=T)
    result <- cbind(count.polygon@data, area=rep(NA, nrow(count.polygon@data)))
    result[,2] <- geosphere::areaPolygon(count.polygon)
# convert from square m to square km
    if (km2 == T) {
        result[,2] <- result[,2]/1000000
        names(result)[2] <- "area.km2"
    }else{
        names(result)[2] <- "area.m2"
    }
    return(result)
}


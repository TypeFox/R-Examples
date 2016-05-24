calcperimeter <- function(x)
{
    ## Verifications
    if (inherits(x, "SpatialPolygonsDataFrame")) {
        d2 <- as(x, "SpatialLinesDataFrame")
    } else {
        d2 <- x
    }
    if (!inherits(d2, "SpatialLinesDataFrame"))
        stop("x should be of class \"SpatialPolygonsDataFrame\" or \"SpatialLinesDataFrame\"")

    uu <- slot(d2, "data")

    di <- unlist(lapply(1:nrow(uu), function(i) {
        coo <- coordinates(d2[i,])
        sum(unlist(lapply(coo, function(y) {
            sum(unlist(lapply(y, function(z) {
                x1 <- z[-1,]
                x2 <- z[-nrow(z),]
                di <- sum(sqrt(((x2[,1]-x1[,1])^2)+((x2[,2]-x1[,2])^2)))
                return(di)
            })))
        })))
    }))

    slot(x,"data") <- data.frame(slot(x,"data"), perimeter=di)
    return(x)
}


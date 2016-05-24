##############################################################################
## package 'secr'
## PG.R
## 2013-11-23
## 2014-08-25 do not require maptools
##############################################################################

PG <- function (CH, poly = NULL, includeNULL = FALSE, plt = FALSE, ...) {
    ## 2014-08-25 require(maptools)
    if (is.null(poly)) {
        poly <- buffer.contour (traps(CH), plt = plt, ...)
        clean <- function(x) {
             ## x may be list or data frame, depending on concave or convex
            if (is.data.frame(x))
                as.matrix(x)
            else  ## assume list
                cbind(x$x, x$y)
        }
        poly <- lapply(poly, clean)
        srl <- lapply(poly, function(x) Polygon(x))
        tmp <- Polygons(srl,1)
        poly <- SpatialPolygons(list(tmp))
        poly <- SpatialPolygonsDataFrame(poly, data = data.frame(ID = names(poly)))
    }
    inpoly <- function (xy) {
        if (is.null(xy)) {
            xy <- matrix(0,nrow=1,ncol=2)
            ip <- NA
        }
        else if (nrow(xy) == 0) {
            ip <- NA
        }
        else {
            ip <- pointsInPolygon(xy, poly)
            if (plt) points (xy[,1], xy[,2], pch=c(1,16)[ip+1])
        }
        sum(ip)/nrow(xy)
    }
    if (detector(traps(CH)) %in% c('polygon','polygonX', 'transect','transectX')) {
        xyl <- split(xy(CH), factor(animalID(CH), levels=rownames(CH) ))  ## keeps order?
    }
    else {
        xyl <- telemetryxy(CH, includeNULL=includeNULL)
    }
    sapply(xyl, inpoly)
}

##############################################################################


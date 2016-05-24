#' @title Percent area sample
#' @description Creates a point sample of polygons where n is based on percent area
#'
#' @param x sp SpatialPolygonsDataFrame object
#' @param pct Percent of area sampled
#' @param join Join polygon attributed to point sample
#' @param msamp Minimum samples
#' @param sf Scaling factor (default is meters to acres conversion factor)
#' @param stype Sampling type ('random', 'regular', 'nonaligned', 'hexagonal')
#' @param ... Additional arguments passed to spsample
#'
#' @note
#' This function results in an adaptive sample based on the area of each polygon
#'
#' @note Depends: sp,  
#'
#' @return A SpatialPointsDataFrame with polygon samples
#'
#' @author Jeffrey S. Evans  <jeffrey_evans@@tnc.org>
#'
#' @examples 
#' require(sp)
#' sr1=Polygons(list(Polygon(cbind(c(180114, 180553, 181127, 181477, 181294, 181007, 180409, 
#'   180162, 180114), c(332349, 332057, 332342, 333250, 333558, 333676, 
#'   332618, 332413, 332349)))),'1')
#' sr2=Polygons(list(Polygon(cbind(c(180042, 180545, 180553, 180314, 179955, 179142, 179437, 
#'   179524, 179979, 180042), c(332373, 332026, 331426, 330889, 330683, 
#'   331133, 331623, 332152, 332357, 332373)))),'2')
#' sr=SpatialPolygons(list(sr1,sr2))
#' srdf=SpatialPolygonsDataFrame(sr, data.frame(row.names=c('1','2'), PIDS=1:2))
#' 
#' ars <- parea.sample(srdf, pct=0.20, stype='random') 
#'   plot(srdf)
#'     plot(ars, pch=20, add=TRUE)
#'
#' @export
parea.sample <- function(x, pct = 0.1, join = FALSE, msamp = 1, sf = 4046.86, stype = "hexagonal", ...) {
    if (!inherits(x, "SpatialPolygonsDataFrame")) 
        stop("Must be a SpatialPolygonsDataFrame object")
    pids <- rownames(x@data)
    psub <- x[rownames(x@data) == pids[1], ]
    ac <- sapply(methods::slot(psub, "polygons"), function(i) methods::slot(i, "area"))/sf
    ns <- round((ac * pct), digits = 0)
    if (ns < msamp) {
        ns <- msamp
    }
    psamp <- sp::spsample(psub, n = ns, type = stype, ...)
    results <- sp::SpatialPointsDataFrame(psamp, data = data.frame(ID = rep(as.numeric(pids[1]), dim(sp::coordinates(psamp))[1])))
    if (length(pids) > 1) {
        for (i in 2:length(pids)) {
            psub <- x[rownames(x@data) == pids[i], ]
            ac <- sapply(methods::slot(psub, "polygons"), function(i) methods::slot(i, "area"))/sf
            ns <- round((ac * pct), digits = 0)
            if (ns < msamp) {
                ns <- msamp
            }
            psamp <- sp::spsample(psub, n = ns, type = stype, ...)
            psamp <- sp::SpatialPointsDataFrame(psamp, data = data.frame(ID = rep(as.numeric(pids[i]), dim(sp::coordinates(psamp))[1])))
            results <- rbind(results, psamp)
        }
    }
    if (join == TRUE) {
        x@data <- cbind(xxID = as.numeric(as.character(rownames(x@data))), x@data)
        results@data <- data.frame(results@data, x@data[match(results@data$ID, x@data$xxID), ])
        x@data <- x@data[, -1]
        return(results)
    } else {
        return(results)
    }
} 

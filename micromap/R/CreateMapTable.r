create_map_table <- function (tmp.map, IDcolumn = NA, poly.thresh=.0001) 
{
    if(is.logical(poly.thresh)){
	 if(!poly.thresh) poly.thresh <- 0
    }

    tot.area <- 0
    lPoly <- vector("list", length(tmp.map@polygons))
    for (i in 1:length(lPoly)) {
        tmp <- tmp.map@polygons[[i]]
        lPoly[[i]] <- lapply(1:length(tmp@Polygons), function(j) cbind(i, 
            j, tmp@Polygons[[j]]@labpt[1], tmp@Polygons[[j]]@labpt[2], 
            tmp@Polygons[[j]]@coords, tmp@Polygons[[j]]@hole, 
            tmp@Polygons[[j]]@area))
	
	  tot.area <- tot.area + tmp@area
    }


    dPoly <- NULL
    for (i in 1:length(lPoly)) {
        for (j in 1:length(lPoly[[i]])) {
            if(lPoly[[i]][[j]][1,8] > tot.area*poly.thresh) dPoly <- rbind(dPoly, lPoly[[i]][[j]])
        }
    }


    if (!is.na(IDcolumn)) {
        tmp.map@data$ID <- tmp.map@data[, IDcolumn]
    }
    else {
        if (is.null(tmp.map@data$ID))  stop("\n ! No ID column exists\n ! No column has been specified as a unique polygon ID\n ! Create.map.table is confused and offended")
    }


    dPoly2 <- data.frame(tmp.map@data$ID[dPoly[, 1]], dPoly)
    names(dPoly2) <- c("ID", "region", "poly", "lab.x", "lab.y", "coordsx", "coordsy", "hole", "area")
    dPoly2 <- transform(dPoly2, poly = (region - 1) * max(dPoly2$poly) + poly)

    tholes <- unique(subset(dPoly2, hole == 1)[, c("poly", "lab.x", "lab.y", "area")])
    tholes$rad <- tholes$area/(2 * pi)
    tnot <- unique(subset(dPoly2, hole == 0)[, c("poly", "lab.x", "lab.y")])


    w1 <- NULL
    if (nrow(tholes) > 0) {
        w <- sapply(1:nrow(tholes), function(j) which.min((tholes[j, 
            2] - tnot[, 2])^2 + (tholes[j, 3] - tnot[, 3])^2))
        w1 <- w[((tholes[, 2] - tnot[w, 2])^2 + (tholes[, 3] - 
            tnot[w, 3])^2) < abs(tholes$rad)]
    }

    dPoly2$plug <- (dPoly2$poly %in% tnot$poly[w1]) * (1 - dPoly2$hole)
    dPoly2 <- dPoly2[, -c(4, 5, 9)]
    dPoly2
}
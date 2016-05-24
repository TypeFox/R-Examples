SpGeoCodH <- function(x, areanames = NULL, occ.thresh = 0) {
    if (class(x) == "spgeoIN") {
        if (class(x$polygons) == "SpatialPolygons") {
            if ("NA" %in% names(x$polygons)) {
                warning("the polygondata contain a polygon named NA. this can cause problems. Please rename")
            }
        }
        if (class(x$polygons) == "SpatialPolygonsDataFrame") {
          if (!areanames %in% names(x$polygons@data)){
            stop(sprintf("column '%s' not found", "areanames"))
          }  
          nam.test <- as.vector(unlist(x$polygons@data[, areanames]))
            if ("NA" %in% nam.test) {
                warning("the polygondata contain a polygon named NA. this can cause problems. Please rename")
            }
        }
        kkk <- .PipSamp(x, columnname = areanames)
        spsum <- .SpSumH(kkk, occ.thresh = occ.thresh)
        
        if (length(spsum) == 0) {
            namco <- c("identifier", names(x$polygons))
            fill <- matrix(0, nrow = length(unique(kkk$identifier)), ncol = length(names(x$polygons)))
            fill <- data.frame(fill)
            spsum <- data.frame(cbind(as.character(unique(kkk$identifier)), fill))
            names(spsum) <- namco
        }
        
        sppol <- .SpPerPolH(spsum)
        
        
        nc <- subset(kkk, is.na(kkk$homepolygon))
        identifier <- x$identifier[as.numeric(rownames(nc))]
        bb <- x$species_coordinates[as.numeric(rownames(nc)), ]
        miss <- data.frame(identifier, bb)
        names(miss) <- c("identifier", "XCOOR", "YCOOR")
        
        out <- list(identifier_in = x$identifier, species_coordinates_in = x$species_coordinates, polygons = x$polygons, sample_table = kkk, 
            spec_table = spsum, polygon_table = sppol, not_classified_samples = miss, coexistence_classified = "NA", areanam = areanames)
        class(out) <- "spgeoOUT"
        return(out)
    } else {
        stop("Function is only defined for class spgeoIN")
    }
} 

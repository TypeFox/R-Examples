RichnessGrid <- function(x, limits = c(-180, 180, -90, 90), reso, type = c("spnum", "abu")) {
    e <- extent(limits)
    xmin <- e[1]
    xmax <- e[2]
    ymin <- e[3]
    ymax <- e[4]
    if (xmin * xmax > 0) {
        cols <- abs(abs(slot(e, "xmax")) - abs(slot(e, "xmin")))
    } else {
        cols <- abs(abs(slot(e, "xmax")) + abs(slot(e, "xmin")))
    }
    if (ymin * ymax > 0) {
        rows <- abs(abs(slot(e, "ymax")) - abs(slot(e, "ymin")))
    } else {
        rows <- abs(abs(slot(e, "ymax")) + abs(slot(e, "ymin")))
    }
    if (class(x) == "spgeoOUT" | class(x) == "spgeoIN") {
        dum <- data.frame(identifier = x$identifier_in, x$species_coordinates_in)
        x <- dum
    }
    if (class(x) == "character" & length(grep(".txt", x)) == 0) {
        if (!requireNamespace("rgbif", quietly = TRUE)) {
          stop("rgbif needed for species name option. Please install it.",
               call. = FALSE)
        } 
        splist <- strsplit(x, " ")
        coords <- rgbif::occ_search(scientificName = x, return = "data", 
                                    limit = 200000, hasCoordinate = T, spatialIssues = F,
                                    fields = c("species", "decimalLongitude","decimalLatitude"))
        coords <- do.call("rbind", coords)
        names(coords) <- c("identifier", "XCOOR", "YCOOR")
        coords <- data.frame(coords[complete.cases(coords),])
        x <- coords
        warning(paste(dim(inp)[1], "geo-referenced records found in GBIF. No data cleaning was performed", sep = " "))
    }
    if (class(x) == "character" & length(grep(".txt", x)) > 0) {
        inp <- read.table(x, sep = "\t", header = T)
        names(inp) <- c("identifier", "XCOOR", "YCOOR")
        x <- inp
    }
    reso <- 60/reso
    ras <- raster(e, ncols = cols * reso, nrows = rows * reso, crs = CRS("+proj=longlat +datum=WGS84"))
    ras <- setValues(ras, 0)
    inp <- split(x, f = x$identifier)
    if (type == "spnum") {
        rast <- lapply(inp, function(x) .rasterSum(x, ras, "div"))
    }
    if (type == "abu") {
        rast <- lapply(inp, function(x) .rasterSum(x, ras, "abu"))
    }
    out <- Reduce("+", rast)
    out[out == 0] <- NA
    
    projection(out) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
    return(out)
}
 

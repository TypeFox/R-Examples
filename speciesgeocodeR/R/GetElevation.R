GetElevation <- function(x) {
    if (class(x) == "data.frame") {
        inp <- x
    }
    if (class(x) == "spgeoIN" | class(x) == "spgeoOUT") {
        inp <- data.frame(identifier = x$identifier, XCOOR = x$species_coordinates[, 1], YCOOR = x$species_coordinates[, 2])
    }
    if (class(x) == "character" & length(grep(".txt", x)) == 0) {
      if (!requireNamespace("rgbif", quietly = TRUE)) {
        stop("rgbif needed for species name option. Please install it.",
             call. = FALSE)
      }  
      coords <- rgbif::occ_search(scientificName = x, return = "data", 
                                  limit = 200000, hasCoordinate = T, spatialIssues = F,
                                  fields = c("species", "decimalLongitude","decimalLatitude"))
      coords <- do.call("rbind", coords)
      names(coords) <- c("identifier", "XCOOR", "YCOOR")
      coords <- data.frame(coords[complete.cases(coords),])
        warning(paste(dim(inp)[1], "geo-referenced records found in GBIF; no data cleaning was performed", sep = " "))
    }
    if (class(x) == "character" & length(grep(".txt", x)) > 0) {
        inp <- read.table(x, sep = "\t", header = T)
        names(inp) <- c("identifier", "XCOOR", "YCOOR")
    }
    
    tt <- list()
    for(i in 1:dim(inp)[1]){
      tt[[i]] <- .getEle(inp[i,])
    }
    
    ele.vector <- suppressWarnings(as.numeric(unlist(tt)))
    
    if (class(x) == "character" & length(grep(".txt", x)) == 0) {
        ele.vector <- cbind(inp, ele.vector)
        return(ele.vector)
    } else {
        return(ele.vector)
    }
} 

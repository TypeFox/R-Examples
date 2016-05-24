ReadPoints <- function(x, y, areanames = NA, verbose = FALSE, cleaning = FALSE) {
  res <- list()
  
  if (class(x) != "character" && class(x) != "data.frame") {
    msg <- sprintf("function not defined for class %s", dQuote(class(x)))
    stop(msg)
  }
  
  if (class(x) == "character" & length(grep(".txt", x)) == 0) {
    if (!requireNamespace("rgbif", quietly = TRUE)) {
      stop("rgbif needed for species name option. Please install it.", 
           call. = FALSE)
    }
    coords <- rgbif::occ_search(scientificName = x, return = "data", limit = 2e+05, 
                                hasCoordinate = T, spatialIssues = F, fields = c("species", "decimalLongitude", 
                                                                                 "decimalLatitude"))
    coords <- do.call("rbind", coords)
    names(coords) <- c("identifier", "XCOOR", "YCOOR")
    coords <- data.frame(coords[complete.cases(coords), ])
    
    if (cleaning == T) {
      cleaner <- GeoClean(coords, isna = TRUE, isnumeric = TRUE, coordinatevalidity = TRUE, 
                          containszero = TRUE, zerozero = TRUE, zerozerothresh = 1, latequallong = TRUE, 
                          GBIFhead = FALSE, countrycentroid = TRUE, contthresh = 0.05, 
                          capitalcoords = TRUE, capthresh = 0.1, countrycheck = FALSE, 
                          referencecountries = speciesgeocodeR::countryref, outp = "cleaned")
      warning(sprintf("%s geo-referenced records found in GBIF; %s  coordinates removed by cleaning", 
                      dim(coords)[1], (dim(coords)[1] - dim(cleaner)[1])))
      cleaner <- cleaner[, 1:3]
      rownames(cleaner) <- NULL
      coords <- cleaner
    } else {
      cleaner <- GeoClean(coords, isna = TRUE, isnumeric = TRUE, coordinatevalidity = TRUE, 
                          containszero = FALSE, zerozero = FALSE, latequallong = FALSE, 
                          GBIFhead = FALSE, countrycentroid = FALSE, capitalcoords = FALSE, 
                          countrycheck = FALSE, outp = "cleaned")
      warning(sprintf("%s geo-referenced records found in GBIF; %s non-numeric coordinates removed", 
                      dim(coords)[1], (dim(coords)[1] - dim(cleaner)[1])))
      cleaner <- cleaner[, 1:3]
      rownames(cleaner) <- NULL
      coords <- cleaner
    }
  }
  
  if (class(x) == "character" & length(grep(".txt", x)) > 0) {
    coords <- read.table(x, sep = "\t", header = T, row.names = NULL)
  }
  
  if (class(x) == "data.frame") {
    coords <- x
    rownames(coords) <- 1:dim(coords)[1]
  }
  
  if (class(y) == "character" | class(y) == "data.frame") {
    if (class(y) == "character" & length(grep(".shp", y)) > 0) {
      poly <- maptools::readShapeSpatial(y)
    } else {
      if (class(y) == "character") {
        polycord <- read.table(y, sep = "\t", header = T)
      }
      if (class(y) == "data.frame") {
        polycord <- y
      }
      if (dim(polycord)[2] != 3) {
        stop("Wrong input format;\ninputfile for polygons must be a tab-delimited text file with three columns")
      }
      if (!is.numeric(polycord[, 2]) || !is.numeric(polycord[, 3])) {
        stop("wrong input format:\nInput polygon coordinates (columns 2 and 3) must be numeric.")
      }
      if (!is.character(polycord[, 1]) && !is.factor(polycord[, 1])) {
        warning("polygon identifier (column 1) should be a string or a factor")
      }
      if (max(polycord[, 2]) > 180) {
        warning(sprintf("check polygon input coordinates; file contains longitude values outside possible range in row: \n                      %s\n Coordinates set to maximum: 180.\n", 
                        rownames(polycord[polycord[, 2] > 180, ])))
        polycord[polycord[, 2] > 180, 2] <- 180
      }
      if (min(polycord[, 2]) < -180) {
        warning(paste("check polygon input coordinates. File contains longitude values outside possible range in row: ", 
                      rownames(polycord[polycord[, 2] < -180, ]), "\n", "Coordinates set to minimum: -180", 
                      sep = ""))
        polycord[polycord[, 2] < -180, ] <- -180
        
      }
      if (max(polycord[, 3]) > 90) {
        warning(paste("check polygon input coordinates. File contains latitude values outside possible range in row:", 
                      rownames(polycord[polycord[, 3] > 90, ]), "\n", "Coordinates set to maximum: 90", 
                      sep = ""))
        polycord[polycord[, 3] > 90, 3] <- 90
      }
      if (min(polycord[, 3]) < -90) {
        warning(paste("check polygon input coordinates. File contains latitude values outside possible range in row:", 
                      rownames(polycord[polycord[, 3] < -90, ]), "\n", "Coordinates set to minimum: -90", 
                      sep = ""))
        polycord[polycord[, 3] < -90, 3] <- -90
      }
      poly <- .Cord2Polygon(polycord)
    }
  }
  
  if (class(y) == "SpatialPolygonsDataFrame" | class(y) == "SpatialPolygons") {
    poly <- y
  }
  
  if (dim(coords)[2] != 3) {
    if (all(c("species", "decimalLatitude", "decimalLongitude") %in% names(coords))) {
      coords <- data.frame(identifier = coords$species, XCOOR = coords$decimalLongitude, 
                           YCOOR = coords$decimalLatitude)
      warning("more than 3 columns in point input. Assuming GBIF file. \nidentifier set to species, XCOOR set to decimalLongitude, YCOOR set to decimalLatitude")
    } else {
      stop(paste("wrong input format: \n", "Inputfile for coordinates must have three columns", 
                 sep = ""))
    }
  }
  
  if (!is.numeric(coords[, 2]) || !is.numeric(coords[, 3])) {
    stop(paste("wrong input format: \n", "Input point coordinates (columns 2 and 3) must be numeric", 
               sep = ""))
  }
  
  if (max(coords[, 2]) > 180) {
    warning(paste("longitude values outside possible range in row:", 
                  rownames(coords[coords[, 2] > 180, ]), ". ", "Row deleted", sep = ""))
    coords <- coords[!coords[, 2] > 180, ]
  }
  if (min(coords[, 2]) < -180) {
    warning(paste("longitude values outside possible range in row: ", 
                  rownames(coords[coords[, 2] < -180, ]), ". ", "Row deleted", sep = ""))
    coords <- coords[!coords[, 2] < -180, ]
  }
  if (max(coords[, 3]) > 90) {
    warning(paste("latitude values outside possible range in row:", 
                  rownames(coords[coords[, 3] > 90, ]), ". ", "Row deleted", sep = ""))
    coords <- coords[!coords[, 3] > 90, ]
  }
  if (min(coords[, 3]) < -90) {
    warning(paste("latitude values outside possible range in row:", 
                  rownames(coords[coords[,3] < -90, ]), ". ", "Row deleted", sep = ""))
    coords <- coords[!coords[, 3] < -90, ]
  }
  if (!is.character(coords[, 1]) && !is.factor(coords[, 1])) {
    warning("coordinate identifier (column 1) should be a string or a factor")
  }
  coords[, 1] <- as.factor(coords[, 1])
  coordi <- coords[, c(2, 3)]
  names(coordi) <- c("XCOOR", "YCOOR")
  
  areanam <- areanames
  
  res <- list(identifier = coords[, 1], species_coordinates = coordi, polygons = poly, 
              areanam = areanam)
  class(res) <- "spgeoIN"
  return(res)
  
} 
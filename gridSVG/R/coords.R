gridSVGCoordsGen <- function() {
  coords <- NULL
  function(newcoords = NULL) {
    if (is.null(newcoords)) {
      coords
    } else if (length(newcoords) == 1 && is.na(newcoords)) {
      coords <<- NULL # Wipe existing info if there is only 1 NA
    } else {
      # Keep existing information, but overwrite any existing definitions
      # associated with any given viewport names
      keepOldName <- ! names(coords) %in% names(newcoords)
      if (any(keepOldName)) {
        existingNames <- names(coords)[keepOldName]
        for (i in 1:length(existingNames)) {
          newcoords[[existingNames[i]]] <- coords[[existingNames[i]]]
        }
      }
      coords <<- newcoords
    }
  }
}

gridSVGCoords <- gridSVGCoordsGen()

validCoordsInfo <- function(vpname) {
  currCoords <- gridSVGCoords()
  if (is.null(currCoords)) {
    stop("No coordinates data has been loaded.")
  } else if (is.null(currCoords[[vpname]])) {
    stop("Viewport not found in coordinates data")
  } else {
    currCoords
  }
}

readCoordsJS <- function(filename) {
  jsData <- readLines(filename)
  jsData <- gsub("var gridSVGCoords = ", "", jsData)
  jsonData <- gsub(";$", "", jsData)
  fromJSON(paste0(jsonData, collapse = "\n"))
}

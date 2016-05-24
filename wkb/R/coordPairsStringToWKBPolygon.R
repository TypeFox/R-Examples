# obj is character string containing delimited coordinate pairs
# by default the delimeter is whitespace
# by default the coordinates in each pair are separated by a comma
# by default reverse is FALSE, meaning that the order of the coordinates is "x,y"
# for geospatial coordinate pairs in the order "latitude,longitude", use reverse = TRUE
# the last coordinate pair is assumed to be identical to the first coordinate pair
coordPairsStringToWKBPolygon <- function(obj, endian = "little", delimiter = "\\s+", separator = ",", reverse = FALSE) {
  wkb <- lapply(X = obj, FUN = function(coords) {
    coords <- strsplit(coords, delimiter)[[1]]
    rc <- rawConnection(raw(0), "r+")
    on.exit(close(rc))
    if(endian == "big") {
      writeBin(as.raw(0L), rc)
    } else {
      writeBin(as.raw(1L), rc)
    }
    writeBin(3L, rc, size = 4, endian = endian)
    writeBin(1L, rc, size = 4, endian = endian)
    writeBin(length(coords), rc, size = 4 , endian = endian)
    lapply(X = coords, FUN = function(coord) {
        coord <- as.numeric(strsplit(coord, separator)[[1]])
        writeBin(coord[1 + reverse], rc, size = 8, endian = endian)
        writeBin(coord[2 - reverse], rc, size = 8, endian = endian)
        NULL
    })
    rawConnectionValue(rc)
  })
  if(identical(version$language, "TERR")) {
    attr(wkb, "SpotfireColumnMetaData") <-
      list(ContentType = "application/x-wkb", MapChart.ColumnTypeId = "Geometry")
  }
  I(wkb)
}

# example usage:
ds <- data.frame(
  name = c(
    "Almost Colorado",
    "Almost Wyoming"
  ),
  geometry = coordPairsStringToWKBPolygon(
    c(
      "41,-109.05 41,-102.04 37,-102.04 37,-109.05 41,-109.05",
      "45,-111.05 45,-104.05 41,-104.05 41,-111.05 45,-111.05"
    ),
    reverse=TRUE
  )
)

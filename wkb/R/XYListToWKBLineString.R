XYListToWKBLineString <- function(obj, endian) {
  wkb <- lapply(obj, function(coords) {
    rc <- rawConnection(raw(0), "r+")
    on.exit(close(rc))
    if(endian == "big") {
      writeBin(as.raw(0L), rc)
    } else {
      writeBin(as.raw(1L), rc)
    }
    writeBin(2L, rc, size = 4, endian = endian)
    writeBin(length(coords[["x"]]), rc, size = 4 , endian = endian)
    mapply(as.numeric(coords[["x"]]), as.numeric(coords[["y"]]), FUN=function(x, y) {
      writeBin(x, rc, size = 8, endian = endian)
      writeBin(y, rc, size = 8, endian = endian)
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

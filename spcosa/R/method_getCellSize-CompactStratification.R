setMethod(
    f = "getCellSize",
    signature = signature(
        object = "CompactStratification"
    ),
    definition = function(object) {
        spatialPixelsDataFrame <- suppressWarnings(as(object, "SpatialPixelsDataFrame"))
        cellSize <- getCellSize(spatialPixelsDataFrame)
        cellSize
    }
)

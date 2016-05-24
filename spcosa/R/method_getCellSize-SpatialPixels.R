setMethod(
    f = "getCellSize",
    signature = signature(
        object = "SpatialPixels"
    ),
    definition = function(object) {
        gridTopology <- getGridTopology(object)
        cellSize <- as(gridTopology, "data.frame")$cellsize
        cellSize
    }
)

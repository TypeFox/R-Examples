setMethod(
    f = "getArea",
    signature = signature(
        object = "CompactStratification"
    ),
    definition = function(object) {
        cellSize <- getCellSize(object)
        cellArea <- prod(cellSize)
        nCells <- getNumberOfCells(object)
        nCells * cellArea        
    }
)

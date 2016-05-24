setMethod(
    f = "stratify",
    signature = signature(
        object = "SpatialPolygons"
    ),
    definition = function(object, nStrata, priorPoints = NULL, maxIterations = 1000L, nTry = 1L,
        nGridCells = 2500L, cellSize, equalArea = FALSE, verbose = getOption("verbose")) {

        # coerce 'object' to an instance of class "SpatialPixels"
        if (missing(cellSize)) {
            object <- spsample(x = object, n = nGridCells, type = "regular")
        } else {
            object <- spsample(x = object, cellsize = cellSize, type = "regular")
        }
        suppressWarnings( # suppress warning when 'grid' has empty column/rows in dimension 2
            gridded(object) <- TRUE
        )

        # stratification
        stratify(object = object, nStrata = nStrata, priorPoints = priorPoints,
            maxIterations = maxIterations, nTry = nTry, equalArea = equalArea, verbose = verbose)
    }
)

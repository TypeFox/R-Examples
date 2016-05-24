setMethod(
    f = "stratify",
    signature = signature(
        object = "SpatialGrid"
    ),
    definition = function(object, nStrata, priorPoints = NULL, maxIterations = 1000L, nTry = 1L,
        equalArea = FALSE, verbose = getOption("verbose")) {
        stratify(object = as(object, "SpatialPixels"),
                 nStrata = nStrata, priorPoints = priorPoints, maxIterations = maxIterations,
                 nTry = nTry, equalArea = equalArea, verbose = verbose)
    }
)

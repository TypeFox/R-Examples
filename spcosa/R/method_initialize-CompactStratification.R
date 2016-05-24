setMethod(
    f = "initialize",
    signature = signature(
        .Object = "CompactStratification"
    ),
    definition = function(.Object, cells, stratumId, centroids, mssd) {
        nCells <- nrow(coordinates(cells))
        if (nCells != length(stratumId)) {
            stop("length of 'stratumId' should match the number of cells", call. = FALSE)
        }
        .Object@cells <- cells
        .Object@stratumId <- stratumId
        .Object@centroids <- centroids
        .Object@mssd <- mssd
        .Object
    }
)

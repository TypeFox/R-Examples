setMethod(
    f = "initialize",
    signature = signature(
        .Object = "CompactStratificationPriorPoints"
    ),
    definition = function(.Object, cells, stratumId, centroids, mssd, priorPoints) {
        .Object <- callNextMethod(.Object, cells = cells, stratumId = stratumId,
            centroids = centroids, mssd = mssd)
        .Object@priorPoints <- priorPoints
        .Object
    }
)

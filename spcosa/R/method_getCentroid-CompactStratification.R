setMethod(
    f = "getCentroid",
    signature = signature(
        object = "CompactStratification"
    ),
    definition = function(object) {
        object@centroids
    }
)

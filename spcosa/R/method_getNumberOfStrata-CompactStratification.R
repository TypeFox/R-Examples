setMethod(
    f = "getNumberOfStrata",
    signature = signature(
        object = "CompactStratification"
    ),
    definition = function(object) {
        length(unique(object@stratumId))
    }
)

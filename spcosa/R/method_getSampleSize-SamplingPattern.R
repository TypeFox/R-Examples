setMethod(
    f = "getSampleSize",
    signature = signature(
        object = "SamplingPattern"
    ),
    definition = function(object) {
        nrow(coordinates(object@sample))
    }
)

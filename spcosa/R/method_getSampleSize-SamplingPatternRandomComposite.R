setMethod(
    f = "getSampleSize",
    signature = signature(
        object = "SamplingPatternRandomComposite"
    ),
    definition = function(object) {
        length(unique(object@composite))
    }
)

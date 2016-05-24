setMethod(
    f = "show",
    signature = signature(
        object = "SamplingPattern"
    ),
    definition = function(object) {
        sampleSize <- getSampleSize(object)
        cat("Object of class", dQuote(class(object)), "\n")
        cat("sample size:", sampleSize, "\n")
    }
)

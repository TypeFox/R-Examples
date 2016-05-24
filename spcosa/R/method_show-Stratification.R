setMethod(
    f = "show",
    signature = signature(
        object = "Stratification"
    ),
    definition = function(object) {
    cat("Object of class", dQuote(class(object)), "\n")
        cat("number of strata:", getNumberOfStrata(object), "\n")
    }
)

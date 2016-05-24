setMethod(
    f = "show",
    signature = signature(
        object = "Statistic"
    ),
    definition = function(object) {
    cat("Object of class", dQuote(class(object)), "\n")
        cat(object@description, "\n")
    }
)

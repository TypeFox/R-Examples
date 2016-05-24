setMethod(
    f = "getObjectiveFunctionValue",
    signature = signature(
        object = "CompactStratification"
    ),
    definition = function(object) {
        object@mssd
    }
)

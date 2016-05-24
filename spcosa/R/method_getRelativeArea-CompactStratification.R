setMethod(
    f = "getRelativeArea",
    signature = signature(
        object = "CompactStratification"
    ),
    definition = function(object) {
        area <- getArea(object)
        area / sum(area)
    }
)

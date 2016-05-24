## Generating function
InfRobRegTypeModel <- function(center = L2RegTypeFamily(), 
                                neighbor = ContNeighborhood()) {
    new("InfRobRegTypeModel", center = center, neighbor = neighbor)
}

## Replace methods
setReplaceMethod("neighbor", "InfRobRegTypeModel", 
    function(object, value){ 
        object@neighbor <- value
        validObject(object) 
        object
    })

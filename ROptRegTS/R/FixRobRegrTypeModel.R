## Generating function
FixRobRegTypeModel <- function(center = RegTypeFamily(), 
                                neighbor = ContNeighborhood()){
    new("FixRobRegTypeModel", center = center, neighbor = neighbor)
}

## Replace methods
setReplaceMethod("neighbor", "FixRobRegTypeModel", 
    function(object, value){ 
        object@neighbor <- value 
        validObject(object)
        object
    })

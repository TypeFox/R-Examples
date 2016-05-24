## access methods
setMethod("name", "Condition", function(object) object@name)

## replace methods
setReplaceMethod("name", "Condition", 
    function(object, value){ object@name <- value; object})

## access methods
setMethod("cond", "UnivariateCondDistribution", function(object) object@cond)

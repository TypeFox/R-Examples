setMethod("name", "RobAStControl", function(object) object@name)
setReplaceMethod("name", "RobAStControl", function(object,value) {object@name <- value; object})

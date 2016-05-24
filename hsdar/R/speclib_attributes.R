setMethod("attribute", signature(object = "Speclib"), 
          function(object)
  return(.attributes(object))
)

setReplaceMethod("attribute", signature(object = "Speclib", value = "matrix"), 
                 function(object, value)
{
  object@attributes <- as.data.frame(value)
  return(object)
}
)

setReplaceMethod("attribute", signature(object = "Speclib", value = "data.frame"),
                 function(object, value)
{
  object@attributes <- value
  return(object)
}
)

.attributes <- function(object)
    return(object@attributes)




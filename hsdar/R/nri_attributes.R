setMethod("attribute", signature(object = "Nri"), 
          function(object)
  return(.attributes(object))
)

setReplaceMethod("attribute", signature(object = "Nri", value = "matrix"), 
                 function(object, value)
{
  object@attributes <- as.data.frame(value)
  return(object)
}
)

setReplaceMethod("attribute", signature(object = "Nri", value = "data.frame"),
                 function(object, value)
{
  object@attributes <- value
  return(object)
}
)




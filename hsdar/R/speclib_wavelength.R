setReplaceMethod("wavelength", signature(object = "Speclib", value = "numeric"), 
                 function(object, value)
{
  object@wavelength <- value
  return(object)
}
)

setMethod("wavelength", signature(object = "Speclib"), 
          function(object)
  return(object@wavelength)
)
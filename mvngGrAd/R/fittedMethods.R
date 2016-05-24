setGeneric("fitted")
setMethod("fitted",
          "movG",
          function(object,...)
          {
            return(object@adjustedPhe)
          }
          )

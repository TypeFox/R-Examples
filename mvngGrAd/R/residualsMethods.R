setGeneric("residuals")
setMethod("residuals",
           "movG",
          function(object,...)
          {
            object <- object@adjModel

            callGeneric(object)
          })

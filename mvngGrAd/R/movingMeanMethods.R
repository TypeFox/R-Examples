setMethod("movingMean",
          signature(object = "movG"),
          function(object,...)
          {
            return(object@movingMean)
          }
          )



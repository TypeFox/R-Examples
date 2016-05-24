setMethod("entryData",
          signature(object = "movG"),
          function(object,...)
          {
            return(round(data.frame(row = object@row,
                              column = object@col,
                              adjustedPhe = object@adjustedPhe,
                              observedPhe = object@observedPhe,
                              movingMean = object@movingMean,
                              nValues = object@nValues),3))
          }
           ) 



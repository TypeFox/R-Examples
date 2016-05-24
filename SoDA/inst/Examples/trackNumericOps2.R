
setMethod("Arith", c("trackNumeric", "numeric"),
          function(e1, e2) {
              xx <- callGeneric(e1@.Data, e2)
              if(length(xx) == length(e1@.Data))
                 e1@.Data <- xx
              else
                 stop("Result of arithmeic operation must preserve original length")
              e1
          })

setMethod("Arith", c("numeric", "trackNumeric"),
          function(e1, e2) {
              xx <- callGeneric(e1, e2@.Data)
              if(length(xx) == length(e2@.Data))
                 e2@.Data <- xx
              else
                 stop("Result of arithmeic operation must preserve original length")
              e2
          })

setMethod("Arith", c("trackNumeric", "trackNumeric"),
          function(e1, e2) {
              if(!(identical(e1@x, e2@x) &&
                 identical(e1@y, e2@y)))
                stop("Objects must have identical x and y coordinates")
              e1@.Data <- callGeneric(e1@.Data, e2@.Data)
              e1
          })


setMethod("Arith", c("trackNumeric", "numeric"),
          function(e1, e2) {
              e1@.Data <- callGeneric(e1@.Data, e2)
              e1
          })

setMethod("Arith", c("numeric", "trackNumeric"),
          function(e1, e2) {
              e2@.Data <- callGeneric(e1, e2@.Data)
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

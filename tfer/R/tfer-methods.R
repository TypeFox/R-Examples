setMethod("show", signature(object = "transfer"),
          function(object)
          print(list(Parameters = para(object), Y = values(object)))
          )


setMethod("[",
          signature(x = "transfer",
                    i = "ANY",
                    j = "missing",
                    drop = "missing"),
          function(x, i, j)
          values(x)[i]
          )


setMethod("Compare",
          signature(e1 = "transfer",
                    e2 = "numeric"),
          function(e1, e2) {
              if (length(e1) < length(e2))
                  stop("incompatible lengths")
              callGeneric(values(e1), as.vector(e2))
          }
          )


setMethod("summary", signature(object = "transfer"),
          function(object){
              return(list(parameters = para(object), values = summary(values(object)),
                          probability = tprob(object)))
          }
          )



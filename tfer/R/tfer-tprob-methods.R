setGeneric("tprob",
           function(object, x)
           standardGeneric("tprob")
           )


setMethod("tprob", signature(object = "transfer",
                             x = "missing"),
          function(object, x){
              return(summary(as.factor(values(object)))/length(values(object)))
          }
          )


setMethod("tprob", signature(object = "transfer",
                             x = "numeric"),
          function(object, x) {
              prob = summary(factor(values(object), levels = 0:max(values(object))))/length(values(object))
              newprob = ifelse(x <= max(values(object)), prob[x+1], 0)
              names(newprob) = x
              return(newprob)
          }
          )

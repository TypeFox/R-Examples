A <- function(object){
  x
}
`A<-` <- function(object,...){
  x
}


# accessor
setMethod(f = "A",
          signature = "IcaSet" ,
          definition = function (object){
              return (object@A)
          }
          )

# set value of attribute A
setReplaceMethod(
      f = "A" ,
      signature = "IcaSet" ,
      definition = function (object, value, valid = TRUE){
        object@A <- value
        if (valid)
            validObject (object)
        return (object)
     }
)

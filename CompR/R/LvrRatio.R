setClass(
  Class="LvrRatio",
  representation=representation(
    Simu="matrix",
    Test="matrix"
    )
)
setGeneric("getSimu",
           function(object)
           {
             standardGeneric("getSimu")
           }
)

setMethod("getSimu","LvrRatio",
          function(object)
          {
            return(object@Simu)
          }
)

setGeneric("getTest",
           function(object)
           {
             standardGeneric("getTest")
           }
)

setMethod("getTest","LvrRatio",
          function(object)
          {
            return(object@Test)
          }
)


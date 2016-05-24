setClass(
  Class="BradleyEstim",
  representation=representation(
    Lvriter="matrix",
    Lvr="numeric",
    Lambda="matrix",
    Pi="list",
    Zh="matrix",
    Ic="matrix",
    Restestglob="list",
    Restestprod="list",
    Varcov="list"
  )
)

setGeneric("getLvriter",
           function(object)
           {
             standardGeneric("getLvriter")
           }
)

setMethod("getLvriter","BradleyEstim",
          function(object)
          {
            return(object@Lvriter)
          }
)

setGeneric("getLvr",
           function(object)
           {
             standardGeneric("getLvr")
           }
)

setMethod("getLvr","BradleyEstim",
          function(object)
          {
            return(object@Lvr)
          }
)

setGeneric("getLambda",
           function(object)
           {
             standardGeneric("getLambda")
           }
)

setMethod("getLambda","BradleyEstim",
          function(object)
          {
            return(object@Lambda)
          }
)

setGeneric("getPi",
           function(object)
           {
             standardGeneric("getPi")
           }
)

setMethod("getPi","BradleyEstim",
          function(object)
          {
            return(object@Pi)
          }
)

setGeneric("getZh",
           function(object)
           {
             standardGeneric("getZh")
           }
)

setMethod("getZh","BradleyEstim",
          function(object)
          {
            return(object@Zh)
          }
)

setGeneric("getIc",
           function(object)
           {
             standardGeneric("getIc")
           }
)

setMethod("getIc","BradleyEstim",
          function(object)
          {
            return(object@Ic)
          }
)

setGeneric("getRestestglob",
           function(object)
           {
             standardGeneric("getRestestglob")
           }
)

setMethod("getRestestglob","BradleyEstim",
          function(object)
          {
            return(object@Restestglob)
          }
)

setGeneric("getRestestprod",
           function(object)
           {
             standardGeneric("getRestestprod")
           }
)

setMethod("getRestestprod","BradleyEstim",
          function(object)
          {
            return(object@Restestprod)
          }
)

setGeneric("getVarcov",
           function(object)
           {
             standardGeneric("getVarcov")
           }
)

setMethod("getVarcov","BradleyEstim",
          function(object)
          {
            return(object@Varcov)
          }
)

setMethod("show","BradleyEstim",
          function(object){
            cat("\n*** Class BradleyEstim, method Show***\n")
            for (i in 1:length(object@Lambda))
            {
              cat("\n******Class",i,"Weight = ",object@Lambda[i]);cat("\n")
              cat("\n***** Bradley's scores******\n")
              print(object@Pi[[i]])
              
            }
          }
)
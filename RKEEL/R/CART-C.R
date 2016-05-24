#Class implementing a Classification Algorithm
  #Implements the CART-C KEEL classification algorithm

CART_C <- function(train, test, maxDepth=90){
  alg <- RKEEL::R6_CART_C$new()
  alg$setParameters(train, test, maxDepth)
  return (alg)
}

R6_CART_C <- R6::R6Class("R6_CART_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #max depth
    maxDepth = 90,


    #Public functions

    #Initialize function
    setParameters = function(train, test, maxDepth=90){

      super$setParameters(train, test)

      self$maxDepth <- maxDepth

    }

  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "Clas-CART.jar",

    #algorithm name
    algorithmName = "CART-C",

    #String with algorithm name
    algorithmString = "Clas-CART",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "maxDepth = ", self$maxDepth, "\n")

      return(text)

    }
  )
)

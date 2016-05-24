#Class implementing a Regression Algorithm
  #Implements the CART-R KEEL regression algorithm

CART_R <- function(train, test, maxDepth=90){
  alg <- RKEEL::R6_CART_R$new()
  alg$setParameters(train, test, maxDepth)
  return (alg)
}

R6_CART_R <- R6::R6Class("R6_CART_R",

  inherit = RegressionAlgorithm,

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
    jarName = "Regr-CART.jar",

    #algorithm name
    algorithmName = "CART-R",

    #String with algorithm name
    algorithmString = "Regr-CART",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "maxDepth = ", self$maxDepth, "\n")

      return(text)

    }
  )
)

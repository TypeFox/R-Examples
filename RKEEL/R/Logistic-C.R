#Class implementing a Classification Algorithm
  #Implements the Logistic-C KEEL classification algorithm

Logistic_C <- function(train, test, ridge=1e-8, maxIter=-1){
  alg <- RKEEL::R6_Logistic_C$new()
  alg$setParameters(train, test, ridge, maxIter)
  return (alg)
}

R6_Logistic_C <- R6::R6Class("R6_Logistic_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #k: num of nearest neighbours
    ridge = 1e-8,

    #distance type
    maxIter = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, ridge=1e-8, maxIter=-1){

      super$setParameters(train, test)

      self$ridge <- ridge
      self$maxIter <- maxIter

    }

  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "Logistic.jar",

    #algorithm name
    algorithmName = "Logistic-C",

    #String with algorithm name
    algorithmString = "Multinomial logistic regression model",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "Ridge = ", self$ridge, "\n")
      text <- paste0(text, "MaxIter = ", self$maxIter, "\n")

      return(text)

    }

  )
)

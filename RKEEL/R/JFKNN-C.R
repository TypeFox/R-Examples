#Class implementing a Classification Algorithm
  #Implements the JFKNN-C KEEL classification algorithm

JFKNN_C <- function(train, test){
  alg <- RKEEL::R6_JFKNN_C$new()
  alg$setParameters(train, test)
  return (alg)
}

R6_JFKNN_C <- R6::R6Class("R6_JFKNN_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties


    #Public functions

    #Initialize function
    setParameters = function(train, test){
      super$setParameters(train, test)

    }

  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "JFKNN.jar",

    #algorithm name
    algorithmName = "JFKNN-C",

    #String with algorithm name
    algorithmString = "Jozwik Fuzzy K Nearest Neighbors Classifier",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""

      return(text)

    }

  )
)

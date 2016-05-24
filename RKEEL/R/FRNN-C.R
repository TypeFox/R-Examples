#Class implementing a Classification Algorithm
  #Implements the FRNN-C KEEL classification algorithm

FRNN_C <- function(train, test){
  alg <- RKEEL::R6_FRNN_C$new()
  alg$setParameters(train, test)
  return (alg)
}

R6_FRNN_C <- R6::R6Class("R6_FRNN_C",

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
    jarName = "FRNN.jar",

    #algorithm name
    algorithmName = "FRNN-C",

    #String with algorithm name
    algorithmString = "Fuzzy Rough Nearest Neighbors Algorithm",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""

      return(text)

    }

  )
)

#Class implementing a Classification Algorithm
  #Implements the CenterNN-C KEEL classification algorithm

CenterNN_C <- function(train, test){
  alg <- RKEEL::R6_CenterNN_C$new()
  alg$setParameters(train, test)
  return (alg)
}

R6_CenterNN_C <- R6::R6Class("R6_CenterNN_C",

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
    jarName = "CenterNN.jar",

    #algorithm name
    algorithmName = "CenterNN-C",

    #String with algorithm name
    algorithmString = "Center based Nearest Neigbor classifier",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""

      return(text)

    }
  )
)

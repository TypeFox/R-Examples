#Class implementing a Classification Algorithm
  #Implements the NM-C KEEL classification algorithm

NM_C <- function(train, test){
  alg <- RKEEL::R6_NM_C$new()
  alg$setParameters(train, test)
  return (alg)
}

R6_NM_C <- R6::R6Class("R6_NM_C",

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
    jarName = "NM.jar",

    #algorithm name
    algorithmName = "NM-C",

    #String with algorithm name
    algorithmString = "Nearest Mean Classifier",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""

      return(text)

    }
  )
)

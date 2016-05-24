#Class implementing a Preprocess Algorithm
  #Implements the MostCommon-MV KEEL preprocess algorithm

MostCommon_MV <- function(train, test){
  alg <- RKEEL::R6_MostCommon_MV$new()
  alg$setParameters(train, test)
  return (alg)
}

R6_MostCommon_MV <- R6::R6Class("R6_MostCommon_MV",

  inherit = PreprocessAlgorithm,

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
    jarName = "MV-MostCommon.jar",

    #algorithm name
    algorithmName = "MostCommon-MV",

    #String with algorithm name
    algorithmString = "Most Common Attribute Value",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""

      return(text)

    }
  )
)

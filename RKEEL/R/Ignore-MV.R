#Class implementing a Preprocess Algorithm
  #Implements the Ignore-MV KEEL preprocess algorithm

Ignore_MV <- function(train, test){
  alg <- RKEEL::R6_Ignore_MV$new()
  alg$setParameters(train, test)
  return (alg)
}

R6_Ignore_MV <- R6::R6Class("R6_Ignore_MV",

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
    jarName = "MV-Ignore.jar",

    #algorithm name
    algorithmName = "Ignore-MV",

    #String with algorithm name
    algorithmString = "Delete Instances with Missing Values",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""

      return(text)

    }
  )
)

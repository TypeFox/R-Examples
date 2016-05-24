#Class implementing a Preprocess Algorithm
  #Implements the CleanAttributes-TR KEEL preprocess algorithm

CleanAttributes_TR <- function(train, test){
  alg <- RKEEL::R6_CleanAttributes_TR$new()
  alg$setParameters(train, test)
  return (alg)
}

R6_CleanAttributes_TR <- R6::R6Class("R6_CleanAttributes_TR",

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
    jarName = "CleanAttributes-TR.jar",

    #algorithm name
    algorithmName = "CleanAttributes-TR",

    #String with algorithm name
    algorithmString = "Cleaning of Redundant or not useful Attributes",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""

      return(text)

    }
  )
)

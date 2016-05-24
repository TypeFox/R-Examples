#Class implementing a Preprocess Algorithm
  #Implements the ZScore-TR KEEL preprocess algorithm

ZScore_TR <- function(train, test){
  alg <- RKEEL::R6_ZScore_TR$new()
  alg$setParameters(train, test)
  return (alg)
}

R6_ZScore_TR <- R6::R6Class("R6_ZScore_TR",

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
    jarName = "Trns-ZScore.jar",

    #algorithm name
    algorithmName = "ZScore-TR",

    #String with algorithm name
    algorithmString = "Nominal To Binary Transformation",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""

      return(text)

    }
  )
)

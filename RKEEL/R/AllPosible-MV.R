#Class implementing a Preprocess Algorithm
  #Implements the AllPosible-MV KEEL preprocess algorithm

AllPosible_MV <- function(train, test){
  alg <- RKEEL::R6_AllPosible_MV$new()
  alg$setParameters(train, test)
  return (alg)
}

R6_AllPosible_MV <- R6::R6Class("R6_AllPosible_MV",

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
    jarName = "MV-AllPossible.jar",

    #algorithm name
    algorithmName = "AllPossible-MV",

    #String with algorithm name
    algorithmString = "Assign all possible values of the attribute",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""

      return(text)

    }
  )
)

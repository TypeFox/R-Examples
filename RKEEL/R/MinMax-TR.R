#Class implementing a Preprocess Algorithm
  #Implements the MinMax-TR KEEL preprocess algorithm

MinMax_TR <- function(train, test, newMin=0.0, newMax=1.0){
  alg <- RKEEL::R6_MinMax_TR$new()
  alg$setParameters(train, test, newMin, newMax)
  return (alg)
}

R6_MinMax_TR <- R6::R6Class("R6_MinMax_TR",

  inherit = PreprocessAlgorithm,

  public = list(

    #Public properties

    #New Min
    newMin = 0.0,

    #New Max
    newMax = 1.0,

    #Public functions

    #Initialize function
    setParameters = function(train, test, newMin=0.0, newMax=1.0){

      super$setParameters(train, test)

      self$newMin <- newMin
      self$newMax <- newMax

    }


  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "Trns-MinMax.jar",

    #algorithm name
    algorithmName = "MinMax-TR",

    #String with algorithm name
    algorithmString = "Min max Normalization",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "newMin = ", self$newMin, "\n")
      text <- paste0(text, "newMax = ", self$newMax, "\n")

      return(text)

    }
  )
)

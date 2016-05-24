#Class implementing a Preprocess Algorithm
  #Implements the DecimalScaling-TR KEEL preprocess algorithm

DecimalScaling_TR <- function(train, test){
  alg <- RKEEL::R6_DecimalScaling_TR$new()
  alg$setParameters(train, test)
  return (alg)
}

R6_DecimalScaling_TR <- R6::R6Class("R6_DecimalScaling_TR",

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
    jarName = "Trns-DecimalScaling.jar",

    #algorithm name
    algorithmName = "DecimalScaling-TR",

    #String with algorithm name
    algorithmString = "Decimal Scaling Normalization",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""

      return(text)

    }
  )
)

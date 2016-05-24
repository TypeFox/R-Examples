#Class implementing a Preprocess Algorithm
  #Implements the Nominal2Binary-TR KEEL preprocess algorithm

Nominal2Binary_TR <- function(train, test){
  alg <- RKEEL::R6_Nominal2Binary_TR$new()
  alg$setParameters(train, test)
  return (alg)
}

R6_Nominal2Binary_TR <- R6::R6Class("R6_Nominal2Binary_TR",

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
    jarName = "Trns-Nominal2Binary.jar",

    #algorithm name
    algorithmName = "Nominal2Binary-TR",

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

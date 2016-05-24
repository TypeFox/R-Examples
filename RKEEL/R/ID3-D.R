#Class implementing a Preprocess Algorithm
  #Implements the ID3-D KEEL preprocess algorithm

ID3_D <- function(train, test){
  alg <- RKEEL::R6_ID3_D$new()
  alg$setParameters(train, test)
  return (alg)
}

R6_ID3_D <- R6::R6Class("R6_ID3_D",

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
    jarName = "Disc-ID3.jar",

    #algorithm name
    algorithmName = "ID3-D",

    #String with algorithm name
    algorithmString = "ID3 Discretizer",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""

      return(text)

    }
  )
)

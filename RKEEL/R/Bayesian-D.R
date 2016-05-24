#Class implementing a Preprocess Algorithm
  #Implements the Bayesian-D KEEL preprocess algorithm

Bayesian_D <- function(train, test){
  alg <- RKEEL::R6_Bayesian_D$new()
  alg$setParameters(train, test)
  return (alg)
}

R6_Bayesian_D <- R6::R6Class("R6_Bayesian_D",

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
    jarName = "Disc-Bayesian.jar",

    #algorithm name
    algorithmName = "Bayesian-D",

    #String with algorithm name
    algorithmString = "Bayesian Discretizer",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""

      return(text)

    }
  )
)

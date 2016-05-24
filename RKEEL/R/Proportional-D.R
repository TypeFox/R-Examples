#Class implementing a Preprocess Algorithm
  #Implements the Proportional-D KEEL preprocess algorithm

Proportional_D <- function(train, test, seed=-1){
  alg <- RKEEL::R6_Proportional_D$new()
  alg$setParameters(train, test, seed)
  return (alg)
}

R6_Proportional_D <- R6::R6Class("R6_Proportional_D",

  inherit = PreprocessAlgorithm,

  public = list(

    #Public properties

    #seed
    seed = -1,

    #Public functions

    #Initialize function
    setParameters = function(train, test, seed=-1){

      super$setParameters(train, test)

      if(seed == -1) {
        self$seed <- sample(1:1000000, 1)
      }
      else {
        self$seed <- seed
      }

    }


  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "Disc-Proportional.jar",

    #algorithm name
    algorithmName = "Proportional-D",

    #String with algorithm name
    algorithmString = "Proportional Discretizer",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")

      return(text)

    }
  )
)

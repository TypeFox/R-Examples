#Class implementing a Classification Algorithm
  #Implements the RBFN-C KEEL classification algorithm

RBFN_C <- function(train, test, neurons=50, seed=-1){
  alg <- RKEEL::R6_RBFN_C$new()
  alg$setParameters(train, test, neurons, seed)
  return (alg)
}

R6_RBFN_C <- R6::R6Class("R6_RBFN_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #num of neurons
    neurons = 50,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, neurons=50, seed=-1){

      super$setParameters(train, test)

      self$neurons <- neurons

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
    jarName = "RBFN_CL.jar",

    #algorithm name
    algorithmName = "RBFN-C",

    #String with algorithm name
    algorithmString = "Radial Basis Function Neural Network for Classification Problems",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "neurons = ", self$neurons, "\n")

      return(text)

    }

  )
)

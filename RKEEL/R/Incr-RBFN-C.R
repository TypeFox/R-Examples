#Class implementing a Classification Algorithm
  #Implements the Incr-RBFN-C KEEL classification algorithm

IncrRBFN_C <- function(train, test, epsilon=0.1, alfa=0.3, delta=0.5, seed=-1){
  alg <- RKEEL::R6_IncrRBFN_C$new()
  alg$setParameters(train, test, epsilon, alfa, delta, seed)
  return (alg)
}

R6_IncrRBFN_C <- R6::R6Class("R6_IncrRBFN_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #epsilon
    epsilon = 0.1,

    #alfa
    alfa = 0.3,

    #delta
    delta = 0.5,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, epsilon=0.1, alfa=0.3, delta=0.5,
                          seed=-1){

      super$setParameters(train, test)

      self$epsilon <- epsilon
      self$alfa <- alfa
      self$delta <- delta

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
    jarName = "RBFN_incremental_CL.jar",

    #algorithm name
    algorithmName = "Incr-RBFN-C",

    #String with algorithm name
    algorithmString = "Incremental Radial Basis Function Neural Network for Classification Problems",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "epsilon = ", self$epsilon, "\n")
      text <- paste0(text, "alfa = ", self$alfa, "\n")
      text <- paste0(text, "delta = ", self$delta, "\n")

      return(text)

    }

  )
)

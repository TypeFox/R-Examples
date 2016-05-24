#Class implementing a Classification Algorithm
  #Implements the FURIA-C KEEL classification algorithm

FURIA_C <- function(train, test, optimizations=2, folds=3, seed=-1){
  alg <- RKEEL::R6_FURIA_C$new()
  alg$setParameters(train, test, optimizations, folds, seed)
  return (alg)
}

R6_FURIA_C <- R6::R6Class("R6_FURIA_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #Number of optimizations
    optimizations = 2,

    #Number of folds
    folds = 3,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, optimizations=2, folds=3, seed=-1){

      super$setParameters(train, test)

      self$optimizations <- optimizations
      self$folds <- folds

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
    jarName = "FURIA-C.jar",

    #algorithm name
    algorithmName = "FURIA-C",

    #String with algorithm name
    algorithmString = "FURIA",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "Number of optimizations = ", self$optimizations, "\n")
      text <- paste0(text, "Number of folds = ", self$folds, "\n")

      return(text)

    }

  )
)

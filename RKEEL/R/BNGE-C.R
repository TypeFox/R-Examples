#Class implementing a Classification Algorithm
  #Implements the BNGE-C KEEL classification algorithm

BNGE_C <- function(train, test, seed=-1){
  alg <- RKEEL::R6_BNGE_C$new()
  alg$setParameters(train, test, seed)
  return (alg)
}

R6_BNGE_C <- R6::R6Class("R6_BNGE_C",

  inherit = ClassificationAlgorithm,

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
    jarName = "BNGE.jar",

    #algorithm name
    algorithmName = "BNGE-C",

    #String with algorithm name
    algorithmString = "BNGE",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")

      return(text)

    }

  )
)

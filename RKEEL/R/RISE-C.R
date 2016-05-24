#Class implementing a Classification Algorithm
  #Implements the RISE-C KEEL classification algorithm

RISE_C <- function(train, test, Q=1, S=2){
  alg <- RKEEL::R6_RISE_C$new()
  alg$setParameters(train, test, Q, S)
  return (alg)
}

R6_RISE_C <- R6::R6Class("R6_RISE_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #Q
    Q = 1,

    #S
    S = 2,


    #Public functions

    #Initialize function
    setParameters = function(train, test, Q=1, S=2){

      super$setParameters(train, test)

      self$Q <- Q
      self$S <- S

    }

  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "RISE.jar",

    #algorithm name
    algorithmName = "RISE-C",

    #String with algorithm name
    algorithmString = "RISE",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "Q = ", self$Q, "\n")
      text <- paste0(text, "S = ", self$S, "\n")

      return(text)

    }

  )
)

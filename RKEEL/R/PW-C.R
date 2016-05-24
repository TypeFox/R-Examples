#Class implementing a Classification Algorithm
  #Implements the PW-C KEEL classification algorithm

PW_C <- function(train, test, beta=8.0, ro=0.001, epsilon=0.001){
  alg <- RKEEL::R6_PW_C$new()
  alg$setParameters(train, test, beta, ro, epsilon)
  return (alg)
}

R6_PW_C <- R6::R6Class("R6_PW_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #Beta
    beta = 8.0,

    #Ro
    ro = 0.001,

    #Epsilon
    epsilon = 0.001,


    #Public functions

    #Initialize function
    setParameters = function(train, test, beta=8.0, ro=0.001, epsilon=0.001){

      super$setParameters(train, test)

      self$beta <- beta
      self$ro <- ro
      self$epsilon <- epsilon

    }


  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "PW.jar",

    #algorithm name
    algorithmName = "PW-C",

    #String with algorithm name
    algorithmString = "Prototipe dependant classifier",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "Beta = ", self$beta, "\n")
      text <- paste0(text, "Ro = ", self$ro, "\n")
      text <- paste0(text, "Epsilon = ", self$epsilon, "\n")

      return(text)

    }
  )
)

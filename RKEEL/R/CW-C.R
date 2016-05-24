#Class implementing a Classification Algorithm
  #Implements the CW-C KEEL classification algorithm

CW_C <- function(train, test, beta=8.0, mu=0.001, epsilon=0.001){
  alg <- RKEEL::R6_CW_C$new()
  alg$setParameters(train, test, beta, mu, epsilon)
  return (alg)
}

R6_CW_C <- R6::R6Class("R6_CW_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #Beta
    beta = 8.0,

    #Mu
    mu = 0.001,

    #Epsilon
    epsilon = 0.001,


    #Public functions

    #Initialize function
    setParameters = function(train, test, beta=8.0, mu=0.001, epsilon=0.001){

      super$setParameters(train, test)

      self$beta <- beta
      self$mu <- mu
      self$epsilon <- epsilon

    }


  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "CW.jar",

    #algorithm name
    algorithmName = "CW-C",

    #String with algorithm name
    algorithmString = "Class dependant classifier",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "Beta = ", self$beta, "\n")
      text <- paste0(text, "Mu = ", self$mu, "\n")
      text <- paste0(text, "Epsilon = ", self$epsilon, "\n")

      return(text)

    }
  )
)

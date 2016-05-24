#Class implementing a Classification Algorithm
  #Implements the CPW-C KEEL classification algorithm

CPW_C <- function(train, test, beta=8.0, mu=0.001, ro=0.001, epsilon=0.001){
  alg <- RKEEL::R6_CPW_C$new()
  alg$setParameters(train, test, beta, mu, ro, epsilon)
  return (alg)
}

R6_CPW_C <- R6::R6Class("R6_CPW_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #Beta
    beta = 8.0,

    #Mu
    mu = 0.001,

    #Ro
    ro = 0.001,

    #Epsilon
    epsilon = 0.001,


    #Public functions

    #Initialize function
    setParameters = function(train, test, beta=8.0, mu=0.001, ro=0.001,
                          epsilon=0.001){

      super$setParameters(train, test)

      self$beta <- beta
      self$mu <- mu
      self$ro <- ro
      self$epsilon <- epsilon

    }


  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "CPW.jar",

    #algorithm name
    algorithmName = "CPW-C",

    #String with algorithm name
    algorithmString = "Class-Prototipe dependant classifier",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "Beta = ", self$beta, "\n")
      text <- paste0(text, "Mu = ", self$mu, "\n")
      text <- paste0(text, "Ro = ", self$ro, "\n")
      text <- paste0(text, "Epsilon = ", self$epsilon, "\n")

      return(text)

    }
  )
)

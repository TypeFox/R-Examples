#Class implementing a Classification Algorithm
  #Implements the Deeps-C KEEL classification algorithm

Deeps_C <- function(train, test, beta=0.12){
  alg <- RKEEL::R6_Deeps_C$new()
  alg$setParameters(train, test, beta)
  return (alg)
}

R6_Deeps_C <- R6::R6Class("R6_Deeps_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #Beta
    beta = 0.12,


    #Public functions

    #Initialize function
    setParameters = function(train, test, beta=0.12){

      super$setParameters(train, test)

      self$beta <- beta

    }


  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "Deeps.jar",

    #algorithm name
    algorithmName = "Deeps-C",

    #String with algorithm name
    algorithmString = "Database extraction of Emerging Patterns classifier",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "Beta = ", self$beta, "\n")

      return(text)

    }
  )
)

#Class implementing a Classification Algorithm
  #Implements the FuzzyNPC-C KEEL classification algorithm

FuzzyNPC_C <- function(train, test, M=2.0){
  alg <- RKEEL::R6_FuzzyNPC_C$new()
  alg$setParameters(train, test, M=2.0)
  return (alg)
}

R6_FuzzyNPC_C <- R6::R6Class("R6_FuzzyNPC_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #M
    M = 2.0,


    #Public functions

    #Initialize function
    setParameters = function(train, test, M=2.0){

      super$setParameters(train, test)

      self$M <- M

    }

  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "FuzzyNPC.jar",

    #algorithm name
    algorithmName = "FuzzyNPC-C",

    #String with algorithm name
    algorithmString = "Fuzzy K Nearest Prototype Classifier",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "M = ", self$M, "\n")

      return(text)

    }

  )
)

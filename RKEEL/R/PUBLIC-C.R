#Class implementing a Classification Algorithm
  #Implements the PUBLIC-C KEEL classification algorithm

PUBLIC_C <- function(train, test, nodesBetweenPrune=25, estimateToPrune="PUBLIC(1)"){
  alg <- RKEEL::R6_PUBLIC_C$new()
  alg$setParameters(train, test, nodesBetweenPrune, estimateToPrune)
  return (alg)
}

R6_PUBLIC_C <- R6::R6Class("R6_PUBLIC_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #nodes between prune
    nodesBetweenPrune = 25,

    #estimate to prune
    estimateToPrune = "PUBLIC(1)",


    #Public functions

    #Initialize function
    setParameters = function(train, test, nodesBetweenPrune=25,
                          estimateToPrune="PUBLIC(1)"){

      super$setParameters(train, test)

      self$nodesBetweenPrune <- nodesBetweenPrune

      if((tolower(estimateToPrune) == "public(1)")  || (tolower(estimateToPrune) == "public(s)") ||
         (tolower(estimateToPrune) == "public(v)")) {
        self$estimateToPrune <- toupper(estimateToPrune)
      }
      else {
        #default value of estimate to prune
        self$estimateToPrune <- "PUBLIC(1)"
      }

    }

  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "PUBLIC.jar",

    #algorithm name
    algorithmName = "PUBLIC-C",

    #String with algorithm name
    algorithmString = "PUBLIC",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "nodesBetweenPrune = ", self$nodesBetweenPrune, "\n")
      text <- paste0(text, "estimateToPrune = ", self$estimateToPrune, "\n")

      return(text)

    }
  )
)

#Class implementing a Classification Algorithm
  #Implements the C45-C KEEL classification algorithm

C45_C <- function(train, test, pruned=TRUE, confidence=0.25, instancesPerLeaf=2){
  alg <- RKEEL::R6_C45_C$new()
  alg$setParameters(train, test, pruned, confidence, instancesPerLeaf)
  return (alg)
}

R6_C45_C <- R6::R6Class("R6_C45_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #pruned
    pruned = TRUE,

    #confidence
    confidence = 0.25,

    #instances per leaf
    instancesPerLeaf = 2,


    #Public functions

    #Initialize function
    setParameters = function(train, test, pruned=TRUE, confidence=0.25,
                          instancesPerLeaf=2){

      super$setParameters(train, test)

      self$pruned <- pruned

      self$confidence <- confidence

      self$instancesPerLeaf <- 2

    }

  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "C45.jar",

    #algorithm name
    algorithmName = "C45-C",

    #String with algorithm name
    algorithmString = "C4.5 Decision Tree",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "pruned = ", self$pruned, "\n")
      text <- paste0(text, "confidence = ", self$confidence, "\n")
      text <- paste0(text, "isntancesPerLeaf = ", self$instancesPerLeaf, "\n")

      return(text)

    }
  )
)



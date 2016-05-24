#Class implementing a Classification Algorithm
  #Implements the PART-C of KEEL

PART_C <- function(train, test, confidence=0.25, itemsetsPerLeaf=2){
  alg <- RKEEL::R6_PART_C$new()
  alg$setParameters(train, test, confidence, itemsetsPerLeaf)
  return (alg)
}

R6_PART_C <- R6::R6Class("R6_PART_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #confidence
    confidence = 0.25,

    #itemsetsPerLeaf
    itemsetsPerLeaf = 2,


    #Public functions

    #Initialize function
    setParameters = function(train, test, confidence=0.25, itemsetsPerLeaf=2){

      super$setParameters(train, test)

      self$confidence <- confidence

      self$itemsetsPerLeaf <- itemsetsPerLeaf

    }

  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "PART.jar",

    #algorithm name
    algorithmName = "PART-C",

    #String with algorithm name
    algorithmString = "PART",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "confidence = ", self$confidence, "\n")
      text <- paste0(text, "itemsetsPerLeaf = ", self$itemsetsPerLeaf, "\n")

      return(text)

    }

  )
)

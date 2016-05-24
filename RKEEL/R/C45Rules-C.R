#Class implementing a Classification Algorithm
  #Implements the C45Rules-C of KEEL

C45Rules_C <- function(train, test, confidence=0.25, itemsetsPerLeaf=2, threshold=10, seed=-1){
  alg <- RKEEL::R6_C45Rules_C$new()
  alg$setParameters(train, test, confidence, itemsetsPerLeaf, threshold, seed)
  return (alg)
}

R6_C45Rules_C <- R6::R6Class("R6_C45Rules_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #confidence
    confidence = 0.25,

    #itemsetsPerLeaf
    itemsetsPerLeaf = 2,

    #threshold
    threshold = 10,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, confidence=0.25, itemsetsPerLeaf=2,
                          threshold=10, seed=-1){

      super$setParameters(train, test)

      self$confidence <- confidence

      self$itemsetsPerLeaf <- itemsetsPerLeaf

      self$threshold <- threshold

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
    jarName = "C45Rules.jar",

    #algorithm name
    algorithmName = "C45Rules-C",

    #String with algorithm name
    algorithmString = "C45Rules",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "confidence = ", self$confidence, "\n")
      text <- paste0(text, "itemsetsPerLeaf = ", self$itemsetsPerLeaf, "\n")
      text <- paste0(text, "threshold = ", self$threshold, "\n")

      return(text)

    }

  )
)

#Class implementing a Classification Algorithm
  #Implements the GFS-AdaBoost-C KEEL classification algorithm

GFS_AdaBoost_C <- function(train, test, numLabels=3, numRules=8, seed=-1){
  alg <- RKEEL::R6_GFS_AdaBoost_C$new()
  alg$setParameters(train, test, numLabels, numRules, seed)
  return (alg)
}

R6_GFS_AdaBoost_C <- R6::R6Class("R6_GFS_AdaBoost_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #Number of labels
    numLabels = 3,

    #Number of rules
    numRules = 8,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, numLabels=3, numRules=8, seed=-1){

      super$setParameters(train, test)

      self$numLabels <- numLabels
      self$numRules <- numRules

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
    jarName = "ClassifierFuzzyAdaBoost.jar",

    #algorithm name
    algorithmName = "GFS-AdaBoost-C",

    #String with algorithm name
    algorithmString = "Fuzzy Rule Learning, AdaBoost Algorithm",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "subAlgorithm = ClasifFuzzyAdaBoost", "\n")
      text <- paste0(text, "dataformat = keel", "\n")
      text <- paste0(text, "numlabels = ", self$numLabels, "\n")
      text <- paste0(text, "numrules = ", self$numRules, "\n")
      text <- paste0(text, "outlabel = FAB", "\n")

      return(text)

    }

  )
)

#Class implementing a Classification Algorithm
  #Implements the GFS-LogitBoost-C KEEL classification algorithm

GFS_LogitBoost_C <- function(train, test, numLabels=3, numRules=25, seed=-1){
  alg <- RKEEL::R6_GFS_LogitBoost_C$new()
  alg$setParameters(train, test, numLabels, numRules, seed)
  return (alg)
}

R6_GFS_LogitBoost_C <- R6::R6Class("R6_GFS_LogitBoost_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #Number of labels
    numLabels = 3,

    #Number of rules
    numRules = 25,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, numLabels=3, numRules=25, seed=-1){

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
    jarName = "ClassifierFuzzyLogitBoost.jar",

    #algorithm name
    algorithmName = "GFS-LogitBoost-C",

    #String with algorithm name
    algorithmString = "Fuzzy Rule Learning, LogitBoost Algorithm",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "subAlgorithm = ClasifFuzzyLogitBoost", "\n")
      text <- paste0(text, "dataformat = keel", "\n")
      text <- paste0(text, "numlabels = ", self$numLabels, "\n")
      text <- paste0(text, "numrules = ", self$numRules, "\n")
      text <- paste0(text, "outlabel = FLB", "\n")

      return(text)

    }

  )
)

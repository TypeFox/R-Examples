#Class implementing a Classification Algorithm
  #Implements the VWFuzzyKNN-C KEEL classification algorithm

VWFuzzyKNN_C <- function(train, test, k=3, init_k=3){
  alg <- RKEEL::R6_VWFuzzyKNN_C$new()
  alg$setParameters(train, test, k, init_k)
  return (alg)
}

R6_VWFuzzyKNN_C <- R6::R6Class("R6_VWFuzzyKNN_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #k
    k = 3,

    #k value (init)
    init_k = 3,

    #Public functions

    #Initialize function
    setParameters = function(train, test, k=3, init_k=3){

      super$setParameters(train, test)

      self$k <- k
      self$init_k <- init_k

    }

  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "VWFuzzyKNN.jar",

    #algorithm name
    algorithmName = "VWFuzzyKNN-C",

    #String with algorithm name
    algorithmString = "Variance Weighted Fuzzy K Nearest Neighbors Classifier",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "K Value = ", self$k, "\n")
      text <- paste0(text, "k Value (init) = ", self$init_k, "\n")

      return(text)

    }

  )
)

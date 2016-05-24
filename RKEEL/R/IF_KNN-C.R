#Class implementing a Classification Algorithm
  #Implements the IF_KNN-C KEEL classification algorithm

IF_KNN_C <- function(train, test, K=3, mA=0.6, vA=0.4, mR=0.3, vR=0.7, k=5){
  alg <- RKEEL::R6_IF_KNN_C$new()
  alg$setParameters(train, test, K, mA, vA, mR, vR, k)
  return (alg)
}

R6_IF_KNN_C <- R6::R6Class("R6_IF_KNN_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #K Value
    K = 3,

    #mA
    mA = 0.6,

    #vA
    vA = 0.4,

    #mR
    mR = 0.3,

    #vR
    vR = 0.7,

    #k Value
    k = 5,


    #Public functions

    #Initialize function
    setParameters = function(train, test, K=3, mA=0.6, vA=0.4, mR=0.3, vR=0.7, k=5){

      super$setParameters(train, test)

      self$K <- K
      self$mA <- mA
      self$vA <- vA
      self$mR <- mR
      self$vR <- vR
      self$k <- k

    }

  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "IF_KNN.jar",

    #algorithm name
    algorithmName = "IF_KNN-C",

    #String with algorithm name
    algorithmString = "Intuitionistic Fuzzy K Nearest Neighbors Classifier ",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "K Value = ", self$K, "\n")
      text <- paste0(text, "mA = ", self$mA, "\n")
      text <- paste0(text, "vA = ", self$vA, "\n")
      text <- paste0(text, "mR = ", self$mR, "\n")
      text <- paste0(text, "vR = ", self$vR, "\n")
      text <- paste0(text, "k Value = ", self$k, "\n")

      return(text)

    }

  )
)

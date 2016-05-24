#Class implementing a Classification Algorithm
  #Implements the CFKNN-C KEEL classification algorithm

CFKNN_C <- function(train, test, k=3, alpha=0.6, seed=-1){
  alg <- RKEEL::R6_CFKNN_C$new()
  alg$setParameters(train, test, k, alpha, seed)
  return (alg)
}

R6_CFKNN_C <- R6::R6Class("R6_CFKNN_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #k
    k = 3,

    #alpha
    alpha = 0.6,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, k=3, alpha=0.6, seed=-1){

      super$setParameters(train, test)

      self$k <- k
      self$alpha <- alpha

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
    jarName = "CFKNN.jar",

    #algorithm name
    algorithmName = "CFKNN-C",

    #String with algorithm name
    algorithmString = "Condensed Fuzzy K Nearest Neighbors Classifier",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "K Value = ", self$k, "\n")
      text <- paste0(text, "Alpha = ", self$alpha, "\n")

      return(text)

    }

  )
)

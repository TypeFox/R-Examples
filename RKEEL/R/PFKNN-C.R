#Class implementing a Classification Algorithm
  #Implements the PFKNN-C KEEL classification algorithm

PFKNN_C <- function(train, test, k=3, seed=-1){
  alg <- RKEEL::R6_PFKNN_C$new()
  alg$setParameters(train, test, k, seed)
  return (alg)
}

R6_PFKNN_C <- R6::R6Class("R6_PFKNN_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #k
    k = 3,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, k=3, seed=-1){

      super$setParameters(train, test)

      self$k <- k

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
    jarName = "PFKNN.jar",

    #algorithm name
    algorithmName = "PFKNN-C",

    #String with algorithm name
    algorithmString = "Pruned Fuzzy K Nearest Neighbors Classifier",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "K Value = ", self$k, "\n")

      return(text)

    }

  )
)

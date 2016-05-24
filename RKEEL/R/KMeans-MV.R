#Class implementing a Preprocess Algorithm
  #Implements the KMeans-MV KEEL preprocess algorithm

KMeans_MV <- function(train, test, k=10, error=100, iterations=100, seed=-1){
  alg <- RKEEL::R6_KMeans_MV$new()
  alg$setParameters(train, test, k, error, iterations, seed)
  return (alg)
}

R6_KMeans_MV <- R6::R6Class("R6_KMeans_MV",

  inherit = PreprocessAlgorithm,

  public = list(

    #Public properties

    #k
    k = 10,

    #error
    error = 100,

    #iterations
    iterations = 100,

    #seed
    seed = -1,

    #Public functions

    #Initialize function
    setParameters = function(train, test, k=10, error=100, iterations=100,
                          seed=-1){

      super$setParameters(train, test)

      self$k <- k
      self$error <- error
      self$iterations <- iterations

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
    jarName = "MV-KMeans.jar",

    #algorithm name
    algorithmName = "KMeans-MV",

    #String with algorithm name
    algorithmString = "K-means Imputation",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "k = ", self$k, "\n")
      text <- paste0(text, "error = ", self$error, "\n")
      text <- paste0(text, "iterations = ", self$iterations, "\n")

      return(text)

    }
  )
)

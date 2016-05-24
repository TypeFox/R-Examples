#Class implementing a Preprocess Algorithm
  #Implements the SSGA-Integer-knn-FS KEEL preprocess algorithm

SSGA_Integer_knn_FS <- function(train, test, paramKNN=1, nEval=5000, pop_size=100, numFeatures=3, seed=-1){
  alg <- RKEEL::R6_SSGA_Integer_knn_FS$new()
  alg$setParameters(train, test, paramKNN, nEval, pop_size, numFeatures, seed)
  return (alg)
}

R6_SSGA_Integer_knn_FS <- R6::R6Class("R6_SSGA_Integer_knn_FS",

  inherit = PreprocessAlgorithm,

  public = list(

    #Public properties

    #paramKNN
    paramKNN = 1,

    #Number of evaluations
    nEval = 5000,

    #Population size
    pop_size = 100,

    #Number of features
    numFeatures = 3,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, paramKNN=1, nEval=5000, pop_size=100,
                          numFeatures=3, seed=-1){

      super$setParameters(train, test)

      self$paramKNN <- paramKNN
      self$pop_size <- pop_size
      self$nEval <- nEval
      self$numFeatures <- numFeatures

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
    jarName = "FS-SSGA-Integer-knn.jar",

    #algorithm name
    algorithmName = "SSGA-Integer-knn-FS",

    #String with algorithm name
    algorithmString = "Steady-state GA with integer coding scheme for wrapper feature selection with k-nn",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "paramKNN = ", self$paramKNN, "\n")
      text <- paste0(text, "nEval = ", self$nEval, "\n")
      text <- paste0(text, "popLength = ", self$pop_size, "\n")
      text <- paste0(text, "numberOfFeatures = ", self$numFeatures, "\n")

      return(text)

    }
  )
)

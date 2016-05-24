#Class implementing a Preprocess Algorithm
  #Implements the UniformFrequency-D KEEL preprocess algorithm

UniformFrequency_D <- function(train, test, numIntervals=10, seed=-1){
  alg <- RKEEL::R6_UniformFrequency_D$new()
  alg$setParameters(train, test, numIntervals, seed)
  return (alg)
}

R6_UniformFrequency_D <- R6::R6Class("R6_UniformFrequency_D",

  inherit = PreprocessAlgorithm,

  public = list(

    #Public properties

    #Num of intervals
    numIntervals = 10,

    #seed
    seed = -1,

    #Public functions

    #Initialize function
    setParameters = function(train, test, numIntervals=10, seed=-1){

      super$setParameters(train, test)

      self$numIntervals <- numIntervals

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
    jarName = "Disc-UniformFrequency.jar",

    #algorithm name
    algorithmName = "UniformFrequency-D",

    #String with algorithm name
    algorithmString = "Uniform Frequency Discretizer",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "numIntervals = ", self$numIntervals, "\n")
      text <- paste0(text, "seed = ", self$seed, "\n")

      return(text)

    }
  )
)

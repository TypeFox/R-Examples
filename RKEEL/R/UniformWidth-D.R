#Class implementing a Preprocess Algorithm
  #Implements the UniformWidth-D KEEL preprocess algorithm

UniformWidth_D <- function(train, test, numIntervals=10){
  alg <- RKEEL::R6_UniformWidth_D$new()
  alg$setParameters(train, test, numIntervals)
  return (alg)
}

R6_UniformWidth_D <- R6::R6Class("R6_UniformWidth_D",

  inherit = PreprocessAlgorithm,

  public = list(

    #Public properties

    #Num of intervals
    numIntervals = 10,

    #Public functions

    #Initialize function
    setParameters = function(train, test, numIntervals=10){

      super$setParameters(train, test)

      self$numIntervals <- numIntervals

    }


  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "Disc-UniformWidth.jar",

    #algorithm name
    algorithmName = "UniformWidth-D",

    #String with algorithm name
    algorithmString = "Uniform Width Discretizer",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "numIntervals = ", self$numIntervals, "\n")

      return(text)

    }
  )
)

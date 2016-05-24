#Class implementing a Preprocess Algorithm
  #Implements the SaturationFilter-F KEEL preprocess algorithm

SaturationFilter_F <- function(train, test, seed=-1){
  alg <- RKEEL::R6_SaturationFilter_F$new()
  alg$setParameters(train, test, seed)
  return (alg)
}

R6_SaturationFilter_F <- R6::R6Class("R6_SaturationFilter_F",

  inherit = PreprocessAlgorithm,

  public = list(

    #Public properties

    #seed
    seed = -1,

    #Public functions

    #Initialize function
    setParameters = function(train, test, seed=-1){

      super$setParameters(train, test)

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
    jarName = "Filter-SaturationFilter.jar",

    #algorithm name
    algorithmName = "SaturationFilter-F",

    #String with algorithm name
    algorithmString = "SaturationFilter",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "noiseSensitivity = 0.75",  "\n")

      return(text)

    }
  )
)

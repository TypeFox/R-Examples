#Class implementing a Classification Algorithm
  #Implements the Slipper-C of KEEL

Slipper_C <- function(train, test, grow_pct=0.66, numBoosting=100, seed=-1){
  alg <- RKEEL::R6_Slipper_C$new()
  alg$setParameters(train, test, grow_pct, numBoosting, seed)
  return (alg)
}

R6_Slipper_C <- R6::R6Class("R6_Slipper_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #grow_pct
    grow_pct = 0.66,

    #numBoosting
    numBoosting = 100,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, grow_pct=0.66, numBoosting=100, seed=-1){

      super$setParameters(train, test)

      self$grow_pct <- grow_pct

      self$numBoosting <- numBoosting

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
    jarName = "Slipper.jar",

    #algorithm name
    algorithmName = "Slipper-C",

    #String with algorithm name
    algorithmString = "Slipper",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "ngrow_pct = ", self$grow_pct, "\n")
      text <- paste0(text, "numBoosting = ", self$numBoosting, "\n")

      return(text)

    }

  )
)

#Class implementing a Classification Algorithm
  #Implements the Ripper-C of KEEL

Ripper_C <- function(train, test, grow_pct=0.66, k=2, seed=-1){
  alg <- RKEEL::R6_Ripper_C$new()
  alg$setParameters(train, test, grow_pct, k, seed)
  return (alg)
}

R6_Ripper_C <- R6::R6Class("R6_Ripper_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #grow_pct
    grow_pct = 0.66,

    #k
    k = 2,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, grow_pct=0.66, k=2, seed=-1){

      super$setParameters(train, test)

      self$grow_pct <- grow_pct

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
    jarName = "Ripper.jar",

    #algorithm name
    algorithmName = "Ripper-C",

    #String with algorithm name
    algorithmString = "Ripper",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "ngrow_pct = ", self$grow_pct, "\n")
      text <- paste0(text, "K = ", self$k, "\n")

      return(text)

    }

  )
)

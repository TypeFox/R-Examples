#Class implementing a Classification Algorithm
  #Implements the PRISM-C of KEEL

PRISM_C <- function(train, test, seed=-1){
  alg <- RKEEL::R6_PRISM_C$new()
  alg$setParameters(train, test, seed)
  return (alg)
}

R6_PRISM_C <- R6::R6Class("R6_PRISM_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, seed=-1){

      super$setParameters(train, test)

      #Check for constraints
      stopText <- ""

      if((hasContinuousData(train)) || (hasContinuousData(test))){
        stopText <- paste0(stopText, "Dataset has continuous data and the algorithm does not accept it.\n")
      }

      if(stopText != ""){
        stop(stopText)
      }

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
    jarName = "Prism.jar",

    #algorithm name
    algorithmName = "PRISM-C",

    #String with algorithm name
    algorithmString = "Prism",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")

      return(text)

    }

  )
)

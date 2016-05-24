#Class implementing a Classification Algorithm
  #Implements the PNN-C KEEL classification algorithm

PNN_C <- function(train, test, seed=-1){
  alg <- RKEEL::R6_PNN_C$new()
  alg$setParameters(train, test, seed)
  return (alg)
}

R6_PNN_C <- R6::R6Class("R6_PNN_C",

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

      if((hasMissingValues(train)) || (hasMissingValues(test))){
        stopText <- paste0(stopText, "Dataset has missing values and the algorithm does not accept it.\n")
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
    jarName = "IG-PNN.jar",

    #algorithm name
    algorithmName = "PNN-C",

    #String with algorithm name
    algorithmString = "PNN",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")

      return(text)

    }

  )
)

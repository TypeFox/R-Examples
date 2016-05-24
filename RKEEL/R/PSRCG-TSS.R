#Class implementing a Preprocess Algorithm
  #Implements the PSRCG-TSS KEEL preprocess algorithm

PSRCG_TSS <- function(train, test, distance="Euclidean"){
  alg <- RKEEL::R6_PSRCG_TSS$new()
  alg$setParameters(train, test, distance)
  return (alg)
}

R6_PSRCG_TSS <- R6::R6Class("R6_PSRCG_TSS",

  inherit = PreprocessAlgorithm,

  public = list(

    #Public properties

    #distance function
    distance = "Euclidean",


    #Public functions

    #Initialize function
    setParameters = function(train, test, distance="Euclidean"){

      super$setParameters(train, test)

      if((tolower(distance) == "euclidean") || (tolower(distance) == "hvdm")){
        self$distance <- distance
      }
      else{
        #Default value of distance function
        self$distance <- "Euclidean"
      }

    }


  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "IS-PSRCG.jar",

    #algorithm name
    algorithmName = "PSRCG-TSS",

    #String with algorithm name
    algorithmString = "Prototipe Selection by Relative Certainty Gain",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "Distance Function = ", self$distance, "\n")

      return(text)

    }
  )
)

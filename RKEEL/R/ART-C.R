#Class implementing a Classification Algorithm
  #Implements the ART-C of KEEL

ART_C <- function(train, test){
  alg <- RKEEL::R6_ART_C$new()
  alg$setParameters(train, test)
  return (alg)
}

R6_ART_C <- R6::R6Class("R6_ART_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties


    #Public functions

    #Initialize function
    setParameters = function(train, test){

		super$setParameters(train, test)

      #Check for constraints
      stopText <- ""

      if((hasContinuousData(train)) || (hasContinuousData(test))){
        stopText <- paste0(stopText, "Dataset has continuous data and the algorithm does not accept it.\n")
      }

      if(stopText != ""){
        stop(stopText)
      }
    }

  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "ART.jar",

    #algorithm name
    algorithmName = "ART-C",

    #String with algorithm name
    algorithmString = "ART",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      return(text)
    }

  )
)

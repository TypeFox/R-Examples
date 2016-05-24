#Class implementing a Preprocess Algorithm
  #Implements the POP-TSS KEEL preprocess algorithm

POP_TSS <- function(train, test){
  alg <- RKEEL::R6_POP_TSS$new()
  alg$setParameters(train, test)
  return (alg)
}

R6_POP_TSS <- R6::R6Class("R6_POP_TSS",

  inherit = PreprocessAlgorithm,

  public = list(

    #Public properties

    #Public functions

    #Initialize function
    setParameters = function(train, test){

      super$setParameters(train, test)

      stopText <- ""

      if(is.null(test)){
        if(hasMissingValues(train)){
          stopText <- paste0(stopText, "Dataset has missing values and the algorithm does not accept it.\n")
        }
      }
      else{
        if((hasMissingValues(train)) || (hasMissingValues(test))){
          stopText <- paste0(stopText, "Dataset has missing values and the algorithm does not accept it.\n")
        }
      }

      if(stopText != ""){
        stop(stopText)
      }

    }


  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "IS-POP.jar",

    #algorithm name
    algorithmName = "POP-TSS",

    #String with algorithm name
    algorithmString = "Patterns by Ordered Projections",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""

      return(text)

    }
  )
)

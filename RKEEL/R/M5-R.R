#Class implementing a Regression Algorithm
  #Implements the M5-R KEEL regression algorithm

M5_R <- function(train, test, type="m", pruningFactor=2, unsmoothed=TRUE){
  alg <- RKEEL::R6_M5_R$new()
  alg$setParameters(train, test, type, pruningFactor, unsmoothed)
  return (alg)
}

R6_M5_R <- R6::R6Class("R6_M5_R",

  inherit = RegressionAlgorithm,

  public = list(

    #Public properties

    #Type
    type = "m",

    #Pruning factor
    pruningFactor = 2,

    #Unsmoothed
    unsmoothed = TRUE,


    #Public functions

    #Initialize function
    setParameters = function(train, test, type="m", pruningFactor=2, unsmoothed=TRUE){

      super$setParameters(train, test)

      if((tolower(type) == "r") || (tolower(type) == "l") || (tolower(type) == "m")){
        self$type <- type
      }
      else{
        #Default value of type
        self$type <- "m"
      }

      self$pruningFactor <- pruningFactor
      self$unsmoothed <- unsmoothed

    }

  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "M5.jar",

    #algorithm name
    algorithmName = "M5-R",

    #String with algorithm name
    algorithmString = "M5",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "type = ", self$type, "\n")
      text <- paste0(text, "pruningFactor = ", self$pruningFactor, "\n")
      text <- paste0(text, "unsmoothed = ", self$unsmoothed, "\n")
      text <- paste0(text, "verbosity = 0", "\n")

      return(text)

    }
  )
)

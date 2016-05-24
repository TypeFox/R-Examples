#Class implementing a Regression Algorithm
  #Implements the M5Rules-R KEEL regression algorithm

M5Rules_R <- function(train, test, pruningFactor=2, heuristic="Coverage"){
  alg <- RKEEL::R6_M5Rules_R$new()
  alg$setParameters(train, test, pruningFactor, heuristic)
  return (alg)
}

R6_M5Rules_R <- R6::R6Class("R6_M5Rules_R",

  inherit = RegressionAlgorithm,

  public = list(

    #Public properties

    #Pruning factor
    pruningFactor = 2,

    #Heuristic
    heuristic = "Coverage",


    #Public functions

    #Initialize function
    setParameters = function(train, test, pruningFactor=2, heuristic="Coverage"){

      super$setParameters(train, test)

      if(tolower(heuristic) == "coverage"){
        self$heuristic <- "Coverage"
      }
      else if((tolower(heuristic) == "rms") || (tolower(heuristic) == "mae") ||
           (tolower(heuristic) == "cc")){
        self$heuristic <- toupper(heuristic)
      }
      else{
        #Default value of heuristic
        self$heuristic <- "Coverage"
      }

      self$pruningFactor <- pruningFactor

    }

  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "M5Rules.jar",

    #algorithm name
    algorithmName = "M5Rules-R",

    #String with algorithm name
    algorithmString = "M5Rules",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "pruningFactor = ", self$pruningFactor, "\n")
      text <- paste0(text, "verbosity = 0", "\n")
      text <- paste0(text, "heuristic = ", self$heuristic, "\n")

      return(text)

    }
  )
)

#Class implementing a Regression Algorithm
#Implements the PolQuadraticLMS-R KEEL regression algorithm

PolQuadraticLMS_R <- function(train, test, seed=-1){
  alg <- RKEEL::R6_PolQuadraticLMS_R$new()
  alg$setParameters(train, test, seed)
  return (alg)
}

R6_PolQuadraticLMS_R <- R6::R6Class("R6_PolQuadraticLMS_R",

  inherit = RegressionAlgorithm,

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
    jarName = "ModelPolQuadraticLMS.jar",

    #algorithm name
    algorithmName = "PolQuadraticLMS-R",

    #String with algorithm name
    algorithmString = "LMS Quadratic Regression",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "subAlgorithm = ModelPolQuadraticLMS", "\n")
      text <- paste0(text, "dataformat = keel", "\n")
      text <- paste0(text, "outlabel = MPQLMS", "\n")

      return(text)

    }

  )
)

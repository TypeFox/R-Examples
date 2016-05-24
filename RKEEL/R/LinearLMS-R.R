#Class implementing a Regression Algorithm
  #Implements the LinearLMS KEEL regression algorithm

LinearLMS_R <- function(train, test, seed=-1){
  alg <- RKEEL::R6_LinearLMS_R$new()
  alg$setParameters(train, test, seed)
  return (alg)
}

R6_LinearLMS_R <- R6::R6Class("R6_LinearLMS_R",

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
    jarName = "ModelLinearLMS.jar",

    #algorithm name
    algorithmName = "LinearLMS-R",

    #String with algorithm name
    algorithmString = "LMS Linear Regression",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "subAlgorithm = ModelLinearLMS", "\n")
      text <- paste0(text, "dataformat = keel", "\n")
      text <- paste0(text, "outlabel = MLLMS", "\n")

      return(text)

    }

  )
)

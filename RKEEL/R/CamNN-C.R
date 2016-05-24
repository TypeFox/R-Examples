#Class implementing a Classification Algorithm
  #Implements the CamNN-C KEEL classification algorithm

CamNN_C <- function(train, test, k=1){
  alg <- RKEEL::R6_CamNN_C$new()
  alg$setParameters(train, test, k)
  return (alg)
}

R6_CamNN_C <- R6::R6Class("R6_CamNN_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #k: num of nearest neighbours
    k = 1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, k=1){

      super$setParameters(train, test)

      self$k <- k

    }


  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "CamNN.jar",

    #algorithm name
    algorithmName = "CamNN-C",

    #String with algorithm name
    algorithmString = "CamNN Classifier",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "K Value = ", self$k, "\n")

      return(text)

    }
  )
)

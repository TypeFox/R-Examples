#Class implementing a Classification Algorithm
  #Implements the KSNN-C KEEL classification algorithm

KSNN_C <- function(train, test, k=1){
  alg <- RKEEL::R6_KSNN_C$new()
  alg$setParameters(train, test, k)
  return (alg)
}

R6_KSNN_C <- R6::R6Class("R6_KSNN_C",

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
    jarName = "KSNN.jar",

    #algorithm name
    algorithmName = "KSNN-C",

    #String with algorithm name
    algorithmString = "K Symmetric Nearest Neighbors Classifier",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "K Value = ", self$k, "\n")

      return(text)

    }
  )
)

#Class implementing a Preprocess Algorithm
  #Implements the KNN-MV KEEL preprocess algorithm

KNN_MV <- function(train, test, k=10){
  alg <- RKEEL::R6_KNN_MV$new()
  alg$setParameters(train, test, k)
  return (alg)
}

R6_KNN_MV <- R6::R6Class("R6_KNN_MV",

  inherit = PreprocessAlgorithm,

  public = list(

    #Public properties

    #k
    k = 10,

    #Public functions

    #Initialize function
    setParameters = function(train, test, k=10){

      super$setParameters(train, test)

      self$k <- k

    }


  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "MV-KNN.jar",

    #algorithm name
    algorithmName = "KNN-MV",

    #String with algorithm name
    algorithmString = "Knn Imputation",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "k = ", self$k, "\n")

      return(text)

    }
  )
)

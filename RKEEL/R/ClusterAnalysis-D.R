#Class implementing a Preprocess Algorithm
  #Implements the ClusterAnalysis-D KEEL preprocess algorithm

ClusterAnalysis_D <- function(train, test){
  alg <- RKEEL::R6_ClusterAnalysis_D$new()
  alg$setParameters(train, test)
  return (alg)
}

R6_ClusterAnalysis_D <- R6::R6Class("R6_ClusterAnalysis_D",

  inherit = PreprocessAlgorithm,

  public = list(

    #Public properties

    #Public functions

    #Initialize function
    setParameters = function(train, test){

      super$setParameters(train, test)

    }


  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "Disc-ClusterAnalysis.jar",

    #algorithm name
    algorithmName = "ClusterAnalysis-D",

    #String with algorithm name
    algorithmString = "Cluster Analysis Discretizer",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""

      return(text)

    }
  )
)

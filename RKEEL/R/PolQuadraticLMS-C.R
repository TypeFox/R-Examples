#Class implementing a Classification Algorithm
  #Implements the PolQuadraticLMS-C KEEL classification algorithm

PolQuadraticLMS_C <- function(train, test, seed=-1){
  alg <- RKEEL::R6_PolQuadraticLMS_C$new()
  alg$setParameters(train, test, seed)
  return (alg)
}

R6_PolQuadraticLMS_C <- R6::R6Class("R6_PolQuadraticLMS_C",

  inherit = ClassificationAlgorithm,

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
    jarName = "ClassifierPolQuadraticLMS.jar",

    #algorithm name
    algorithmName = "PolQuadraticLMS-C",

    #String with algorithm name
    algorithmString = "Least Mean Square Quadratic classifier",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "subAlgorithm = ClasifPolQuadraticLMS", "\n")
      text <- paste0(text, "dataformat = keel", "\n")
      text <- paste0(text, "outlabel = CPQLMS", "\n")

      return(text)

    }

  )
)

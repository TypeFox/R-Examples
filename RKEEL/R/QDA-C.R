#Class implementing a Classification Algorithm
  #Implements the QDA-C KEEL classification algorithm

QDA_C <- function(train, test, seed=-1){
  alg <- RKEEL::R6_QDA_C$new()
  alg$setParameters(train, test, seed)
  return (alg)
}

R6_QDA_C <- R6::R6Class("R6_QDA_C",

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
    jarName = "ClassifierADQuadratic.jar",

    #algorithm name
    algorithmName = "QDA-C",

    #String with algorithm name
    algorithmString = "Quadratic Discriminant Analysis",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "subAlgorithm = ClasifADQuadratic", "\n")
      text <- paste0(text, "dataformat = keel", "\n")
      text <- paste0(text, "outlabel = ADQ", "\n")

      return(text)

    }

  )
)

#Class implementing a Preprocess Algorithm
  #Implements the ANR-F KEEL preprocess algorithm

ANR_F <- function(train, test, seed=-1){
  alg <- RKEEL::R6_ANR_F$new()
  alg$setParameters(train, test, seed)
  return (alg)
}

R6_ANR_F <- R6::R6Class("R6_ANR_F",

  inherit = PreprocessAlgorithm,

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
    jarName = "Filter-ANR.jar",

    #algorithm name
    algorithmName = "ANR-F",

    #String with algorithm name
    algorithmString = "ANR",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "subAlgorithm = ClasifMLPerceptron", "\n")
      text <- paste0(text, "dataformat = keel", "\n")
      text <- paste0(text, "outlabel = CMLP", "\n")
      text <- paste0(text, "topologymlp = 1", "\n")

      return(text)

    }
  )
)

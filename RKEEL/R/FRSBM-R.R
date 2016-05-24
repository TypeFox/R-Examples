#Class implementing a Regression Algorithm
  #Implements the FRSBM-R KEEL regression algorithm

FRSBM_R <- function(train, test, numrules=1, sigma=0.0001, seed=-1){
  alg <- RKEEL::R6_FRSBM_R$new()
  alg$setParameters(train, test, numrules, sigma, seed)
  return (alg)
}

FRNN_C <- function(train, test){
  alg <- RKEEL::R6_FRNN_C$new()
  alg$setParameters(train, test)
  return (alg)
}

R6_FRSBM_R <- R6::R6Class("R6_FRSBM_R",

  inherit = RegressionAlgorithm,

  public = list(

    #Public properties

    #number of rules
    numrules = 1,

    #sigma fss98
    sigma = 0.0001,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, numrules=1, sigma=0.0001, seed=-1){

      super$setParameters(train, test)

      self$numrules <- numrules
      self$sigma <- sigma

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
    jarName = "ModelFSS98.jar",

    #algorithm name
    algorithmName = "FRSBM-R",

    #String with algorithm name
    algorithmString = "Fuzzy and Random Sets based Modeling",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "subAlgorithm = ModelFSS98", "\n")
      text <- paste0(text, "dataformat = keel", "\n")
      text <- paste0(text, "outlabel = FSS98", "\n")
      text <- paste0(text, "numrules = ", self$numrules, "\n")
      text <- paste0(text, "sigmafss98 = ", self$sigma, "\n")

      return(text)

    }
  )
)

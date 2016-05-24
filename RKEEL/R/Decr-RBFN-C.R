#Class implementing a Classification Algorithm
  #Implements the Decr-RBFN-R KEEL classification algorithm

DecrRBFN_C <- function(train, test, percent=0.1, num_neurons_ini=20, alfa=0.3, seed=-1){
  alg <- RKEEL::R6_DecrRBFN_C$new()
  alg$setParameters(train, test, percent, num_neurons_ini, alfa, seed)
  return (alg)
}

R6_DecrRBFN_C <- R6::R6Class("R6_DecrRBFN_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #percent
    percent = 0.1,

    #num of initial neurons
    num_neurons_ini = 20,

    #alfa
    alfa = 0.3,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, percent=0.1, num_neurons_ini=20,
                          alfa=0.3, seed=-1){

      super$setParameters(train, test)

      self$percent <- percent
      self$num_neurons_ini <- num_neurons_ini
      self$alfa <- alfa

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
    jarName = "RBFN_decremental_CL.jar",

    #algorithm name
    algorithmName = "Decr-RBFN-C",

    #String with algorithm name
    algorithmString = "Decremental Radial Basis Function Neural Network for classification problems",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "percent = ", self$percent, "\n")
      text <- paste0(text, "nNeuronsIni = ", self$num_neurons_ini, "\n")
      text <- paste0(text, "alfa = ", self$alfa, "\n")

      return(text)

    }

  )
)

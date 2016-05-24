#Class implementing a Regression Algorithm
#Implements the RBFN-R KEEL regression algorithm

RBFN_R <- function(train, test, neurons=50, seed=-1){
  alg <- RKEEL::R6_RBFN_R$new()
  alg$setParameters(train, test, neurons, seed)
  return (alg)
}

R6_RBFN_R <- R6::R6Class("RBFN_R",

  inherit = RegressionAlgorithm,

  public = list(

    #Public properties

    #num of neurons
    neurons = 50,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, neurons=50, seed=-1){

      super$setParameters(train, test)

      self$neurons <- neurons

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
    jarName = "RBFN.jar",

    #algorithm name
    algorithmName = "RBFN-R",

    #String with algorithm name
    algorithmString = "Radial Basis Function Neural Network for Regression Problems",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "neurons = ", self$neurons, "\n")

      return(text)

    }

  )
)

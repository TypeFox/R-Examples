#Class implementing a Classification Algorithm
  #Implements the Kernel-C KEEL classification algorithm

Kernel_C <- function(train, test, sigma=0.01, seed=-1){
  alg <- RKEEL::R6_Kernel_C$new()
  alg$setParameters(train, test, sigma, seed)
  return (alg)
}

R6_Kernel_C <- R6::R6Class("R6_Kernel_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #sigma kernel
    sigma = 0.01,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, sigma=0.01, seed=-1){

      super$setParameters(train, test)

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
    jarName = "ClassifierKernel.jar",

    #algorithm name
    algorithmName = "Kernel-C",

    #String with algorithm name
    algorithmString = "Kernel Classifier",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "subAlgorithm = ClasifKernel", "\n")
      text <- paste0(text, "dataformat = keel", "\n")
      text <- paste0(text, "outlabel = CK", "\n")
      text <- paste0(text, "sigmaKernel = ", self$sigma, "\n")

      return(text)

    }

  )
)

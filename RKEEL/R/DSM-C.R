#Class implementing a Classification Algorithm
  #Implements the DSM-C KEEL classification algorithm

DSM_C <- function(train, test, iterations=100, percentage=10, alpha_0=0.1, seed=-1){
  alg <- RKEEL::R6_DSM_C$new()
  alg$setParameters(train, test, iterations, percentage, alpha_0, seed)
  return (alg)
}

R6_DSM_C <- R6::R6Class("R6_DSM_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #Number of iterations
    iterations = 100,

    #Percentage respect training size
    percentage = 10,

    #alpha_0
    alpha_0 = 0.1,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, iterations=100, percentage=10,
                          alpha_0=0.1, seed=-1){

		  super$setParameters(train, test)

      #Check for constraints
      stopText <- ""

      if((hasMissingValues(train)) || (hasMissingValues(test))){
        stopText <- paste0(stopText, "Dataset has missing values and the algorithm does not accept it.\n")
      }

      if(stopText != ""){
        stop(stopText)
      }

      self$iterations <- iterations
      self$percentage <- percentage
      self$alpha_0 <- alpha_0

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
    jarName = "IG-DSM.jar",

    #algorithm name
    algorithmName = "DSM-C",

    #String with algorithm name
    algorithmString = "DSM",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "number_of_iterations = ", self$iterations, "\n")
      text <- paste0(text, "percentage_respect_training_size = ", self$percentage, "\n")
      text <- paste0(text, "alpha_0 = ", self$alpha_0, "\n")

      return(text)

    }

  )
)

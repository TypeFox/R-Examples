#Class implementing a Classification Algorithm
  #Implements the Bojarczuk_GP-C KEEL classification algorithm

Bojarczuk_GP_C <- function(train, test, population_size=200, max_generations=200, max_deriv_size=20, rec_prob=0.8, copy_prob=0.01, seed=-1){
  alg <- RKEEL::R6_Bojarczuk_GP_C$new()
  alg$setParameters(train, test, population_size, max_generations, max_deriv_size, rec_prob, copy_prob, seed)
  return (alg)
}

R6_Bojarczuk_GP_C <- R6::R6Class("R6_Bojarczuk_GP_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #population size
    population_size = 200,

    #max num of generations
    max_generations = 200,

    #max deriv size
    max_deriv_size = 20,

    #recombination probability
    rec_prob = 0.8,

    #copy probability
    copy_prob = 0.01,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, population_size=200, max_generations=200,
                          max_deriv_size=20, rec_prob=0.8, copy_prob=0.01,
                          seed=-1){

		  super$setParameters(train, test)

      #Check for constraints
      stopText <- ""

      if((hasMissingValues(train)) || (hasMissingValues(test))){
        stopText <- paste0(stopText, "Dataset has missing values and the algorithm does not accept it.\n")
      }

      if(stopText != ""){
        stop(stopText)
      }

      self$population_size <- population_size
      self$max_generations <- max_generations
      self$max_deriv_size <- max_deriv_size
      self$rec_prob <- rec_prob
      self$copy_prob <- copy_prob

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
    jarName = "Bojarczuk_GP.jar",

    #algorithm name
    algorithmName = "Bojarczuk_GP-C",

    #String with algorithm name
    algorithmString = "Bojarczuk_GP",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "population-size = ", self$population_size, "\n")
      text <- paste0(text, "max-generations = ", self$max_generations, "\n")
      text <- paste0(text, "max-deriv-size = ", self$max_deriv_size, "\n")
      text <- paste0(text, "rec-prob = ", self$rec_prob, "\n")
      text <- paste0(text, "copy-prob = ", self$copy_prob, "\n")

      return(text)

    }

  )
)

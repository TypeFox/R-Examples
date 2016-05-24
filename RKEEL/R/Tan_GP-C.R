#Class implementing a Classification Algorithm
  #Implements the Tan_GP-C KEEL classification algorithm

Tan_GP_C <- function(train, test, population_size=150, max_generations=100, max_deriv_size=20, rec_prob=0.8, mut_prob=0.1, copy_prob=0.01, w1=0.7, w2=0.8, elitist_prob=0.06, support=0.03, seed=-1){
  alg <- RKEEL::R6_Tan_GP_C$new()
  alg$setParameters(train, test, population_size, max_generations, max_deriv_size, rec_prob, mut_prob, copy_prob, w1, w2, elitist_prob, support, seed)
  return (alg)
}

R6_Tan_GP_C <- R6::R6Class("Tan_GP_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #population size
    population_size = 150,

    #max num of generations
    max_generations = 100,

    #max deriv size
    max_deriv_size = 20,

    #recombination probability
    rec_prob = 0.8,

    #mutation probability
    mut_prob = 0.1,

    #copy probability
    copy_prob = 0.01,

    #w1
    w1 = 0.7,

    #w2
    w2 = 0.8,

    #elitist probability
    elitist_prob = 0.06,

    #support
    support = 0.03,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, population_size=150, max_generations=100,
                          max_deriv_size=20, rec_prob=0.8, mut_prob=0.1,
                          copy_prob=0.01, w1=0.7, w2=0.8, elitist_prob=0.06,
                          support=0.03, seed=-1){

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
      self$mut_prob <- mut_prob
      self$copy_prob <- copy_prob
      self$w1 <- w1
      self$w2 <- w2
      self$elitist_prob <- elitist_prob
      self$support <- support

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
    jarName = "Tan_GP.jar",

    #algorithm name
    algorithmName = "Tan_GP-C",

    #String with algorithm name
    algorithmString = "Tan_GP",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "population-size = ", self$population_size, "\n")
      text <- paste0(text, "max-generations = ", self$max_generations, "\n")
      text <- paste0(text, "max-deriv-size = ", self$max_deriv_size, "\n")
      text <- paste0(text, "rec-prob = ", self$rec_prob, "\n")
      text <- paste0(text, "mut-prob = ", self$mut_prob, "\n")
      text <- paste0(text, "copy-prob = ", self$copy_prob, "\n")
      text <- paste0(text, "w1 = ", self$w1, "\n")
      text <- paste0(text, "w2 = ", self$w2, "\n")
      text <- paste0(text, "elitist-prob = ", self$elitist_prob, "\n")
      text <- paste0(text, "support = ", self$support, "\n")

      return(text)

    }

  )
)

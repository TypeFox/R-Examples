#Class implementing a Classification Algorithm
  #Implements the FCRA-C of KEEL

FCRA_C <- function(train, test, generations=50, pop_size=30, length_S_C=10, WCAR=10.0, WV=1.0, crossover_prob=1.0, mut_prob=0.01, n1=0.001, n2=0.1, max_iter=100, linguistic_values=5, seed=-1){
  alg <- RKEEL::R6_FCRA_C$new()
  alg$setParameters(train, test, generations, pop_size, length_S_C, WCAR, WV, crossover_prob, mut_prob, n1, n2, max_iter, linguistic_values, seed)
  return (alg)
}

R6_FCRA_C <- R6::R6Class("R6_FCRA_C",

  inherit = AssociativeClassificationAlgorithm,

  public = list(

    #Public properties

    #Number of generations
    generations = 50,

    #Population size
    pop_size = 30,

    #Length of S and C
    length_S_C = 10,

    #Weight of the classification accuracy rate (WCAR)
    WCAR = 10.0,

    #Weight of the number of fuzzy rules (WV)
    WV = 1.0,

    #Crossover probability
    crossover_prob = 1.0,

    #Mutation Probability (per gen)
    mut_prob = 0.01,

    #Learning rate n1 (Nozaki method)
    n1 = 0.001,

    #Learning rate n2 (Nozaki method)
    n2 = 0.1,

    #Number of iterations Jmax (Nozaki method)
    max_iter = 100,

    #Number of Linguistic Values
    linguistic_values = 5,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, generations=50, pop_size=30,
                          length_S_C=10, WCAR=10.0, WV=1.0, crossover_prob=1.0,
                          mut_prob=0.01, n1=0.001, n2=0.1, max_iter=100,
                          linguistic_values=5, seed=-1){

		  super$setParameters(train, test)

      #Check for constraints
      stopText <- ""

      if((hasMissingValues(train)) || (hasMissingValues(test))){
        stopText <- paste0(stopText, "Dataset has missing values and the algorithm does not accept it.\n")
      }

      if(stopText != ""){
        stop(stopText)
      }

      self$generations <- generations
      self$pop_size <- pop_size
      self$length_S_C <- length_S_C
      self$WCAR <- WCAR
      self$WV <- WV
      self$crossover_prob <- crossover_prob
      self$mut_prob <- mut_prob
      self$n1 <- n1
      self$n2 <- n2
      self$max_iter <- max_iter
      self$linguistic_values <- linguistic_values

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
    jarName = "Clas-Fuzzy-FCRA.jar",

    #algorithm name
    algorithmName = "FCRA-C",

    #String with algorithm name
    algorithmString = "Finding fuzzy classification rules based on the Apriori algorithm",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "Number of Generations = ", self$generations, "\n")
      text <- paste0(text, "Population size = ", self$pop_size, "\n")
      text <- paste0(text, "Length of S and C = ", self$length_S_C, "\n")
      text <- paste0(text, "Weight of the classification accuracy rate (WCAR) = ", self$WCAR, "\n")
      text <- paste0(text, "Weight of the number of fuzzy rules (WV) = ", self$WV, "\n")
      text <- paste0(text, "Crossover Probability = ", self$crossover_prob, "\n")
      text <- paste0(text, "Mutation Probability (per gen) = ", self$mut_prob, "\n")
      text <- paste0(text, "Learning rate n1 (Nozaki method) = ", self$n1, "\n")
      text <- paste0(text, "Learning rate n2 (Nozaki method) = ", self$n2, "\n")
      text <- paste0(text, "Number of iterations Jmax (Nozaki method) = ", self$max_iter, "\n")
      text <- paste0(text, "Number of Linguistic Values = ", self$linguistic_values, "\n")

      return(text)

    }

  )
)

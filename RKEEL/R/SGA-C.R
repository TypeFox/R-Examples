#Class implementing a Classification Algorithm
  #Implements the SGA-C KEEL classification algorithm

SGA_C <- function(train, test, mut_prob_1to0=0.01, mut_prob_0to1=0.001, cross_prob=1, pop_size=50, evaluations=10000, alfa=0.5, selection_type="orden_based", k=1, distance="Euclidean", seed=-1){
  alg <- RKEEL::R6_SGA_C$new()
  alg$setParameters(train, test, mut_prob_1to0, mut_prob_0to1, cross_prob, pop_size, evaluations, alfa, selection_type, k, distance, seed)
  return (alg)
}

R6_SGA_C <- R6::R6Class("R6_SGA_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #1 to 0 Mutation Probability = 0.01
    mut_prob_1to0 = 0.01,

    #0 to 1 Mutation Probability = 0.001
    mut_prob_0to1 = 0.001,

    #Cross probability
    cross_prob = 1,

    #Population size
    pop_size = 50,

    #Number of evaluations
    evaluations = 10000,

    #Alfa Equilibrate Factor
    alfa = 0.5,

    #Selection type
    selection_type = "orden_based",

    #Number of Neighbours
    k = 1,

    #Distance Function
    distance = "Euclidean",

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, mut_prob_1to0=0.01, mut_prob_0to1=0.001,
                          cross_prob=1, pop_size=50, evaluations=10000,
                          alfa=0.5, selection_type="orden_based", k=1,
                          distance="Euclidean", seed=-1){

      super$setParameters(train, test)

      self$mut_prob_1to0 <- mut_prob_1to0
      self$mut_prob_0to1 <- mut_prob_0to1
      self$cross_prob <- cross_prob
      self$pop_size <- pop_size
      self$evaluations <- evaluations
      self$alfa <- alfa

      if((tolower(selection_type) == "orden_based") || (tolower(selection_type) == "binary_tournament")){
        self$selection_type <- tolower(selection_type)
      }
      else{
        #Default value of selection type
        self$selection_type <- "orden_based"
      }

      self$k <- k

      if(tolower(distance) == "euclidean"){
        self$distance <- "Euclidean"
      }
      else if(tolower(distance) == "hvdm"){
        self$distance <- "HVDM"
      }
      else{
        #Default value of distance type
        self$distance <- "Euclidean"
      }

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
    jarName = "IS-SGA.jar",

    #algorithm name
    algorithmName = "SGA",

    #String with algorithm name
    algorithmString = "Steady-State Genetic Algorithm for Instance Selection",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "1 to 0 Mutation Probability = ", self$mut_prob_1to0, "\n")
      text <- paste0(text, "0 to 1 Mutation Probability = ", self$mut_prob_0to1, "\n")
      text <- paste0(text, "Cross Probability = ", self$cross_prob, "\n")
      text <- paste0(text, "Population Size = ", self$pop_size, "\n")
      text <- paste0(text, "Number of Evaluations = ", self$evaluations, "\n")
      text <- paste0(text, "Alfa Equilibrate Factor = ", self$alfa, "\n")
      text <- paste0(text, "Selection Type = ", self$selection_type, "\n")
      text <- paste0(text, "Number of Neighbours = ", self$k, "\n")
      text <- paste0(text, "Distance Function = ", self$distance, "\n")

      return(text)

    }

  )
)

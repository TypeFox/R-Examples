#Class implementing a Classification Algorithm
  #Implements the CHC-C KEEL classification algorithm

CHC_C <- function(train, test, pop_size=50, evaluations=10000, alfa=0.5, restart_change=0.35, prob_restart=0.25, prob_diverge=0.05, k=1, distance="Euclidean", seed=-1){
  alg <- RKEEL::R6_CHC_C$new()
  alg$setParameters(train, test, pop_size, evaluations, alfa, restart_change, prob_restart, prob_diverge, k, distance, seed)
  return (alg)
}

R6_CHC_C <- R6::R6Class("R6_CHC_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #Population size
    pop_size = 50,

    #Number of evaluations
    evaluations = 10000,

    #Alfa Equilibrate Factor
    alfa = 0.5,

    #Percentage of Change in Restart
    restart_change = 0.35,

    #0 to 1 Probability in Restart
    prob_restart = 0.25,

    #0 to 1 Probability in Diverge
    prob_diverge = 0.05,

    #Number of Neighbours
    k = 1,

    #Distance Function
    distance = "Euclidean",

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, pop_size=50, evaluations=10000, alfa=0.5,
                          restart_change=0.35, prob_restart=0.25,
                          prob_diverge=0.05, k=1, distance="Euclidean", seed=-1){

      super$setParameters(train, test)

      self$pop_size <- pop_size
      self$evaluations <- evaluations
      self$alfa <- alfa
      self$restart_change <- restart_change
      self$prob_restart <- prob_restart
      self$prob_diverge <- prob_diverge
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
    jarName = "IS-CHC.jar",

    #algorithm name
    algorithmName = "CHC",

    #String with algorithm name
    algorithmString = "CHC Adaptative Search for Instance Selection",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "Population Size = ", self$pop_size, "\n")
      text <- paste0(text, "Number of Evaluations = ", self$evaluations, "\n")
      text <- paste0(text, "Alfa Equilibrate Factor = ", self$alfa, "\n")
      text <- paste0(text, "Percentage of Change in Restart = ", self$restart_change, "\n")
      text <- paste0(text, "0 to 1 Probability in Restart = ", self$prob_restart, "\n")
      text <- paste0(text, "0 to 1 Probability in Diverge = ", self$prob_diverge, "\n")
      text <- paste0(text, "Number of Neighbours = ", self$k, "\n")
      text <- paste0(text, "Distance Function = ", self$distance, "\n")

      return(text)

    }

  )
)

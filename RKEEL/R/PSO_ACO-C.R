#Class implementing a Classification Algorithm
  #Implements the PSO_ACO-C KEEL classification algorithm

PSO_ACO_C <- function(train, test, max_uncovered_samples=20, min_saples_by_rule=2, max_iterations_without_converge=100, enviromentSize=3, numParticles=100, x=0.72984, c1=2.05, c2=2.05, seed=-1){
  alg <- RKEEL::R6_PSO_ACO_C$new()
  alg$setParameters(train, test, max_uncovered_samples, min_saples_by_rule, max_iterations_without_converge, enviromentSize, numParticles, x, c1, c2, seed)
  return (alg)
}

R6_PSO_ACO_C <- R6::R6Class("R6_PSO_ACO_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #Max number of uncovered samples
    max_uncovered_samples = 20,

    #Min number of samples by rule
    min_saples_by_rule = 2,

    #Max iterations without converge
    max_iterations_without_converge = 100,

    #Enviroment size
    enviromentSize = 3,

    #Number of particles
    numParticles = 100,

    #X
    x = 0.72984,

    #C1
    c1 = 2.05,

    #C2
    c2 = 2.05,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test,
                          max_uncovered_samples=20, min_saples_by_rule=2,
                          max_iterations_without_converge=100,
                          enviromentSize=3, numParticles=100, x=0.72984,
                          c1=2.05, c2=2.05, seed=-1){

      super$setParameters(train, test)

      self$max_uncovered_samples <- max_uncovered_samples
      self$min_saples_by_rule <- min_saples_by_rule
      self$max_iterations_without_converge <- max_iterations_without_converge
      self$enviromentSize <- enviromentSize
      self$numParticles <- numParticles
      self$x <- x
      self$c1 <- c1
      self$c2 <- c2

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
    jarName = "PSO_ACO.jar",

    #algorithm name
    algorithmName = "PSO_ACO-C",

    #String with algorithm name
    algorithmString = "PSO/ACO",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "Max_number_of_uncovered_samples = ", self$max_uncovered_samples, "\n")
      text <- paste0(text, "Min_number_of_samples_by_rule = ", self$min_saples_by_rule, "\n")
      text <- paste0(text, "Max_Iterations_without_converge = ", self$max_iterations_without_converge, "\n")
      text <- paste0(text, "Enviroment_Size = ", self$enviromentSize, "\n")
      text <- paste0(text, "Number_of_particles = ", self$numParticles, "\n")
      text <- paste0(text, "X = ", self$x, "\n")
      text <- paste0(text, "C1 = ", self$c1, "\n")
      text <- paste0(text, "C2 = ", self$c2, "\n")
      text <- paste0(text, "flag = 10", "\n")

      return(text)

    }

  )
)

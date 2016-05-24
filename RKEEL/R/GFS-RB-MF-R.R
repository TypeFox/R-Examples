#Class implementing a Regression Algorithm
  #Implements the GFS-RB-MF-R KEEL regression algorithm

GFS_RB_MF_R <- function(train, test, numLabels=3, popSize=50, generations=100, crossProb=0.9, mutProb=0.1, seed=-1){
  alg <- RKEEL::R6_GFS_RB_MF_R$new()
  alg$setParameters(train, test, numLabels, popSize, generations, crossProb, mutProb, seed)
  return (alg)
}

R6_GFS_RB_MF_R <- R6::R6Class("R6_GFS_RB_MF_R",

  inherit = RegressionAlgorithm,

  public = list(

    #Public properties

    #Number of labels
    numLabels = 3,

    #Population size
    popSize = 50,

    #Number of generations
    generations = 100,

    #Crossover probability
    crossProb = 0.9,

    #Mutation probability
    mutProb = 0.1,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, numLabels=3, popSize=50, generations=100,
                          crossProb=0.9, mutProb=0.1, seed=-1){

      super$setParameters(train, test)

      self$numLabels <- numLabels
      self$popSize <- popSize
      self$generations <- generations
      self$crossProb <- crossProb
      self$mutProb <- mutProb

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
    jarName = "GFS_RB_MF.jar",

    #algorithm name
    algorithmName = "GFS-RB-MF-R",

    #String with algorithm name
    algorithmString = "Genetic-Based Fuzzy Rule Base Construction and Membership Functions Tuning",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "Number of Labels = ", self$numLabels, "\n")
      text <- paste0(text, "Population Size = ", self$popSize, "\n")
      text <- paste0(text, "Number of Generations = ", self$generations, "\n")
      text <- paste0(text, "Crossover probability = ", self$crossProb, "\n")
      text <- paste0(text, "Mutation Probability = ", self$mutProb, "\n")

      return(text)

    },

    #Create config files function
    writeKeelConfig = function(){

      #Ficheros de entrada
      inputDataString = paste0("\"../datasets/", private$dataName, "/", private$trainFilename, "\" \"../datasets/", private$dataName, "/", private$trainFilename, "\" \"../datasets/", private$dataName, "/", private$testFilename, "\"")
      #Ficheros de salida
      #dataName <- strsplit(trainData, "/")[[1]][1]
      outputDataString = paste0("\"../results/", private$algorithmName, ".",private$dataName, "/result0", ".tra\" \"../results/", private$algorithmName, ".", private$dataName, "/result0", ".tst\" \"../results/", private$algorithmName, ".", private$dataName, "/result0", "e0.txt\" \"../results/", private$algorithmName, ".", private$dataName, "/result0", "e1.txt\"")

      text <- ""
      text <- paste0(text, "algorithm = ", private$algorithmString)
      text <- paste0(text, "\ninputData = ", inputDataString)
      text <- paste0(text, "\noutputData = ", outputDataString)
      text <- paste0(text, "\n\n")
      text <- paste0(text, private$getParametersText())
      text <- paste0(text, "\n")

      writeLines(text, con = paste0(private$mainPath, "/scripts/", private$algorithmName, "/", private$dataName, "/", "config0s0", ".txt"))

    }
  )
)

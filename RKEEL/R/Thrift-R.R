#Class implementing a Regression Algorithm
  #Implements the Thrift-R KEEL regression algorithm

Thrift_R <- function(train, test, numLabels=3, popSize=61, evaluations=10000, crossProb=0.6, mutProb=0.1, seed=-1){
  alg <- RKEEL::R6_Thrift_R$new()
  alg$setParameters(train, test, numLabels, popSize, evaluations, crossProb, mutProb, seed)
  return (alg)
}

R6_Thrift_R <- R6::R6Class("R6_Thrift_R",

  inherit = RegressionAlgorithm,

  public = list(

    #Public properties

    #Number of labels
    numLabels = 3,

    #Population size
    popSize = 61,

    #Number of evaluations
    evaluations = 10000,

    #Crossover probability
    crossProb = 0.6,

    #Mutation probability
    mutProb = 0.1,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, numLabels=3, popSize=61,
                          evaluations=10000, crossProb=0.6, mutProb=0.1,
                          seed=-1){

      super$setParameters(train, test)

      self$numLabels <- numLabels
      self$popSize <- popSize
      self$evaluations <- evaluations
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
    jarName = "Thrift.jar",

    #algorithm name
    algorithmName = "Thrift-R",

    #String with algorithm name
    algorithmString = "Fuzzy Rule Learning, Thrift Algorithm",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "Number of labels = ", self$numLabels, "\n")
      text <- paste0(text, "Population Size = ", self$popSize, "\n")
      text <- paste0(text, "Number of Evaluations = ", self$evaluations, "\n")
      text <- paste0(text, "Crossover Probability = ", self$crossProb, "\n")
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

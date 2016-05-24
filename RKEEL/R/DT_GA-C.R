#Class implementing a Classification Algorithm
  #Implements the Hybrid Decision Tree/Genetic Algorithm KEEL classification algorithm

DT_GA_C <- function(train, test, confidence=0.25, instancesPerLeaf=2, geneticAlgorithmApproach="GA-LARGE-SN", threshold=10, numGenerations=50, popSize=200, crossoverProb=0.8, mutProb=0.01, seed=-1){
  alg <- RKEEL::R6_DT_GA_C$new()
  alg$setParameters(train, test, confidence, instancesPerLeaf, geneticAlgorithmApproach, threshold, numGenerations, popSize, crossoverProb, mutProb, seed)
  return (alg)
}

R6_DT_GA_C <- R6::R6Class("R6_DT_GA_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #confidence
    confidence = 0.25,

    #instances per leaf
    instancesPerLeaf = 2,

    #Genetic Algorithm Approach
    geneticAlgorithmApproach = "GA-LARGE-SN",

    #Threshold S to consider a Small Disjunt
    threshold = 10,

    #Number of Total Generations for the GA
    numGenerations = 50,

    #Number of chromosomes in the population
    popSize = 200,

    #Crossover Probability
    crossoverProb = 0.8,

    #Mutation Probability
    mutProb = 0.01,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, confidence=0.25, instancesPerLeaf=2,
                          geneticAlgorithmApproach="GA-LARGE-SN", threshold=10,
                          numGenerations=50, popSize=200, crossoverProb=0.8,
                          mutProb=0.01, seed=-1){

		  super$setParameters(train, test)

      #Check for constraints
      stopText <- ""

      if((hasMissingValues(train)) || (hasMissingValues(test))){
        stopText <- paste0(stopText, "Dataset has missing values and the algorithm does not accept it.\n")
      }

      if(stopText != ""){
        stop(stopText)
      }

      self$confidence <- confidence

      self$instancesPerLeaf <- 2

      if((tolower(geneticAlgorithmApproach) == "ga-small") || (tolower(geneticAlgorithmApproach) == "ga-large-sn")){
        self$geneticAlgorithmApproach <- toupper(geneticAlgorithmApproach)
      }
      else{
        #Default value of genetic algorithm approach
        self$geneticAlgorithmApproach <- "GA-LARGE-SN"
      }

      self$threshold <- threshold

      self$numGenerations <- numGenerations

      self$popSize <- popSize

      self$crossoverProb <- crossoverProb

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
    jarName = "DT_GA.jar",

    #algorithm name
    algorithmName = "DT_GA-C",

    #String with algorithm name
    algorithmString = "Hybrid Decision Tree/Genetic Algorithm",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "confidence = ", self$confidence, "\n")
      text <- paste0(text, "isntancesPerLeaf = ", self$instancesPerLeaf, "\n")
      text <- paste0(text, "Genetic Algorithm Approach = ", self$geneticAlgorithmApproach, "\n")
      text <- paste0(text, "Threshold S to consider a Small Disjunt = ", self$threshold, "\n")
      text <- paste0(text, "Number of Total Generations for the GA = ", self$numGenerations, "\n")
      text <- paste0(text, "Number of chromosomes in the population = ", self$popSize, "\n")
      text <- paste0(text, "Crossover Probability = ", self$crossoverProb, "\n")
      text <- paste0(text, "Mutation Probability = ", self$mutProb, "\n")

      return(text)

    }
  )
)

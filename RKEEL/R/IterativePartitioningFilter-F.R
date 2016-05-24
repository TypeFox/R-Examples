#Class implementing a Preprocess Algorithm
  #Implements the IterativePartitioningFilter-F KEEL preprocess algorithm

IterativePartitioningFilter_F <- function(train, test, numPartitions=5, filterType="consensus", confidence=0.25, itemsetsPerLeaf=2, seed=-1){
  alg <- RKEEL::R6_IterativePartitioningFilter_F$new()
  alg$setParameters(train, test, numPartitions, filterType, confidence, itemsetsPerLeaf, seed)
  return (alg)
}

R6_IterativePartitioningFilter_F <- R6::R6Class("R6_IterativePartitioningFilter_F",

  inherit = PreprocessAlgorithm,

  public = list(

    #Public properties

    #Number of partitions
    numPartitions = 5,

    #Filter type
    filterType = "consensus",

    #Confidence
    confidence = 0.25,

    #Itemsets per leaf
    itemsetsPerLeaf = 2,

    #seed
    seed = -1,

    #Public functions

    #Initialize function
    setParameters = function(train, test, numPartitions=5, filterType="consensus",
                          confidence=0.25, itemsetsPerLeaf=2, seed=-1){

      super$setParameters(train, test)

      self$numPartitions <- numPartitions

      if((tolower(filterType) == "consensus") || (tolower(filterType) == "majority")){
        self$filterType <- filterType
      }
      else{
        #Default value of filter type
        self$filterType <- "consensus"
      }

      self$confidence <- confidence
      self$itemsetsPerLeaf <- itemsetsPerLeaf

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
    jarName = "Filter-IterativePartitioningFilter.jar",

    #algorithm name
    algorithmName = "IterativePartitioningFilter-F",

    #String with algorithm name
    algorithmString = "IterativePartitioningFilter",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "numPartitions = ", self$numPartitions, "\n")
      text <- paste0(text, "filterType = ", self$filterType, "\n")
      text <- paste0(text, "prune = true",  "\n")
      text <- paste0(text, "confidence = ", self$confidence, "\n")
      text <- paste0(text, "itemsetsPerLeaf = ", self$itemsetsPerLeaf, "\n")

      return(text)

    }
  )
)

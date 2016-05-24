#Class implementing a Preprocess Algorithm
  #Implements the Relief-FS KEEL preprocess algorithm

Relief_FS <- function(train, test, paramKNN=1, relevanceThreshold=0.20, numInstancesSampled=1000, seed=-1){
  alg <- RKEEL::R6_Relief_FS$new()
  alg$setParameters(train, test, paramKNN, relevanceThreshold, numInstancesSampled, seed)
  return (alg)
}

R6_Relief_FS <- R6::R6Class("R6_Relief_FS",

  inherit = PreprocessAlgorithm,

  public = list(

    #Public properties

    #paramKNN
    paramKNN = 1,

    #Relevance threshold
    relevanceThreshold = 0.20,

    #Num Instances Sampled
    numInstancesSampled = 1000,

    #seed
    seed = -1,

    #Public functions

    #Initialize function
    setParameters = function(train, test, paramKNN=1, relevanceThreshold=0.20,
                          numInstancesSampled=1000, seed=-1){

      super$setParameters(train, test)

      stopText <- ""

      if(is.null(test)){
        if(hasMissingValues(train)){
          stopText <- paste0(stopText, "Dataset has missing values and the algorithm does not accept it.\n")
        }
        if(hasContinuousData(train)){
          stopText <- paste0(stopText, "Dataset has continuous data and the algorithm does not accept it.\n")
        }
      }
      else{
        if((hasMissingValues(train)) || (hasMissingValues(test))){
          stopText <- paste0(stopText, "Dataset has missing values and the algorithm does not accept it.\n")
        }
        if((hasContinuousData(train)) || (hasContinuousData(test))){
          stopText <- paste0(stopText, "Dataset has continuous data and the algorithm does not accept it.\n")
        }
      }

      if(stopText != ""){
        stop(stopText)
      }

      self$paramKNN <- paramKNN
      self$relevanceThreshold <- relevanceThreshold
      self$numInstancesSampled <- numInstancesSampled

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
    jarName = "FS-Relief.jar",

    #algorithm name
    algorithmName = "Relief-FS",

    #String with algorithm name
    algorithmString = "Relief",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "paramKNN = ", self$paramKNN, "\n")
      text <- paste0(text, "relevanceThreshold = ", self$relevanceThreshold, "\n")
      text <- paste0(text, "numInstancesSampled = ", self$numInstancesSampled, "\n")

      return(text)

    }
  )
)

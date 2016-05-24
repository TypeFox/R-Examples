#Class implementing a Classification Algorithm
  #Implements the CFAR-C of KEEL

CFAR_C <- function(train, test, min_support=0.1, min_confidence=0.85, threshold=0.15, num_labels=5, seed=-1){
  alg <- RKEEL::R6_CFAR_C$new()
  alg$setParameters(train, test, min_support, min_confidence, threshold, num_labels, seed)
  return (alg)
}

R6_CFAR_C <- R6::R6Class("R6_CFAR_C",

  inherit = AssociativeClassificationAlgorithm,

  public = list(

    #Public properties

    #Minimum Support
    min_support = 0.1,

    #Minimum Confidence
    min_confidence = 0.85,

    #Cut Threshold (MS)
    threshold = 0.15,

    #Number of labels
    num_labels = 5,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, min_support=0.1, min_confidence=0.85,
                          threshold=0.15, num_labels=5, seed=-1){

		  super$setParameters(train, test)

      #Check for constraints
      stopText <- ""

      if((hasMissingValues(train)) || (hasMissingValues(test))){
        stopText <- paste0(stopText, "Dataset has missing values and the algorithm does not accept it.\n")
      }

      if(stopText != ""){
        stop(stopText)
      }

		  self$min_support <- min_support
      self$min_confidence <- min_confidence
      self$threshold <- threshold
      self$num_labels <- num_labels

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
    jarName = "Clas-Fuzzy-CFAR.jar",

    #algorithm name
    algorithmName = "CFAR-C",

    #String with algorithm name
    algorithmString = "Classification with Fuzzy Association Rules (CFAR)",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "Minimum Support = ", self$min_support, "\n")
      text <- paste0(text, "Minimum Confidence = ", self$min_confidence, "\n")
      text <- paste0(text, "Cut Threshold (MS) = ", self$threshold, "\n")
      text <- paste0(text, "Number of Labels = ", self$num_labels, "\n")

      return(text)

    }

  )
)

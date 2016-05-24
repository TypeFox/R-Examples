#Class implementing a Classification Algorithm
  #Implements the Fuzzy-FARCHD-C of KEEL

FuzzyFARCHD_C <- function(train, test, linguistic_values=5, min_support=0.05, max_confidence=0.8, depth_max=3, K=2, max_evaluations=15000, pop_size=50, alpha=0.15, bits_per_gen=30, inference_type=1, seed=-1){
  alg <- RKEEL::R6_FuzzyFARCHD_C$new()
  alg$setParameters(train, test, linguistic_values, min_support, max_confidence, depth_max, K, max_evaluations, pop_size, alpha, bits_per_gen, inference_type, seed)
  return (alg)
}

R6_FuzzyFARCHD_C <- R6::R6Class("R6_FuzzyFARCHD_C",

  inherit = AssociativeClassificationAlgorithm,

  public = list(

    #Public properties

    #Number of Linguistic Values
    linguistic_values = 5,

    #Minimum Support = 0.05
    min_support = 0.05,

    #Maximum Confidence = 0.8
    max_confidence = 0.8,

    #Depth of the trees (Depthmax) = 3
    depth_max = 3,

    #Parameter K of the prescreening = 2
    K = 2,

    #Maximum number of evaluacions = 15000
    max_evaluations = 15000,

    #Population size = 50
    pop_size = 50,

    #Parameter alpha = 0.15
    alpha = 0.15,

    #Bits per gen = 30
    bits_per_gen = 30,

    #Type of inference = 1
    inference_type = 1,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, linguistic_values=5, min_support=0.05,
                          max_confidence=0.8, depth_max=3, K=2,
                          max_evaluations=15000, pop_size=50, alpha=0.15,
                          bits_per_gen=30, inference_type=1, seed=-1){

		  super$setParameters(train, test)

      #Check for constraints
      stopText <- ""

      if((hasMissingValues(train)) || (hasMissingValues(test))){
        stopText <- paste0(stopText, "Dataset has missing values and the algorithm does not accept it.\n")
      }

      if(stopText != ""){
        stop(stopText)
      }

      self$linguistic_values <- linguistic_values
      self$min_support <- min_support
      self$max_confidence <- max_confidence
      self$depth_max <- depth_max
      self$K <- K
      self$max_evaluations <- max_evaluations
      self$pop_size <- pop_size
      self$alpha <- alpha
      self$bits_per_gen <- bits_per_gen
      self$inference_type <- inference_type

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
    jarName = "Clas-Fuzzy-FARCHD.jar",

    #algorithm name
    algorithmName = "Fuzzy-FARCHD-C",

    #String with algorithm name
    algorithmString = "FARCHD",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "Number of Linguistic Values = ", self$linguistic_values, "\n")
      text <- paste0(text, "Minimum Support = ", self$min_support, "\n")
      text <- paste0(text, "Maximum Confidence = ", self$max_confidence, "\n")
      text <- paste0(text, "Depth of the trees (Depthmax) = ", self$depth_max, "\n")
      text <- paste0(text, "Parameter K of the prescreening = ", self$K, "\n")
      text <- paste0(text, "Maximum number of evaluacions = ", self$max_evaluations, "\n")
      text <- paste0(text, "Population size = ", self$pop_size, "\n")
      text <- paste0(text, "Parameter alpha = ", self$alpha, "\n")
      text <- paste0(text, "Bits per gen = ", self$bits_per_gen, "\n")
      text <- paste0(text, "Type of inference = ", self$inference_type, "\n")

      return(text)

    }

  )
)

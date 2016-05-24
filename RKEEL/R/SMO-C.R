#Class implementing a Classification Algorithm
  #Implements the SMO-C KEEL classification algorithm

SMO_C <- function(train, test, C=1.0, toleranceParameter=0.001, epsilon=1.0e-12, RBFKernel_gamma=0.01, normalized_PolyKernel_exponent=1, normalized_PolyKernel_useLowerOrder=FALSE, PukKernel_omega=1.0, PukKernel_sigma=1.0, StringKernel_lambda=0.5, StringKernel_subsequenceLength=3, StringKernel_maxSubsequenceLength=9, StringKernel_normalize=FALSE,  StringKernel_pruning="None", KernelType="PolyKernel", FitLogisticModels=FALSE, ConvertNominalAttributesToBinary=TRUE, PreprocessType="Normalize", seed=-1){
  alg <- RKEEL::R6_SMO_C$new()
  alg$setParameters(train, test, C, toleranceParameter, epsilon, RBFKernel_gamma, normalized_PolyKernel_exponent, normalized_PolyKernel_useLowerOrder, PukKernel_omega, PukKernel_sigma, StringKernel_lambda, StringKernel_subsequenceLength, StringKernel_maxSubsequenceLength, StringKernel_normalize,  StringKernel_pruning, KernelType, FitLogisticModels, ConvertNominalAttributesToBinary, PreprocessType, seed)
  return (alg)
}

R6_SMO_C <- R6::R6Class("R6_SMO_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #C
    C = 1.0,

    #tolerance parameter
    toleranceParameter = 0.001,

    #epsilon
    epsilon = 1.0e-12,

    #RBFKernel gamma
    RBFKernel_gamma = 0.01,

    #-Normalized PolyKernel exponent
    normalized_PolyKernel_exponent = 1,

    #-Normalized PolyKernel useLowerOrder
    normalized_PolyKernel_useLowerOrder = FALSE,

    #PukKernel Omega
    PukKernel_omega = 1.0,

    #PukKernel sigma
    PukKernel_sigma = 1.0,

    #String Kernel lambda
    StringKernel_lambda = 0.5,

    #String Kernel subsequence length
    StringKernel_subsequenceLength = 3,

    #String Kernel max subsequence length
    StringKernel_maxSubsequenceLength = 9,

    #String Kernel normalize
    StringKernel_normalize = FALSE,

    #String Kernel pruning
    StringKernel_pruning = "None",

    #Kernel type
    KernelType = "PolyKernel",

    #Fit Logistic Models
    FitLogisticModels = FALSE,

    #Convert Nominal Attributes To Binary
    ConvertNominalAttributesToBinary = TRUE,

    #Preprocess type
    PreprocessType = "Normalize",

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test,
                          C=1.0, toleranceParameter=0.001, epsilon=1.0e-12,
                          RBFKernel_gamma=0.01, normalized_PolyKernel_exponent=1,
                          normalized_PolyKernel_useLowerOrder=FALSE,
                          PukKernel_omega=1.0, PukKernel_sigma=1.0,
                          StringKernel_lambda=0.5,
                          StringKernel_subsequenceLength=3,
                          StringKernel_maxSubsequenceLength=9,
                          StringKernel_normalize=FALSE,
                          StringKernel_pruning="None", KernelType="PolyKernel",
                          FitLogisticModels=FALSE,
                          ConvertNominalAttributesToBinary=TRUE,
                          PreprocessType="Normalize", seed=-1){

		  super$setParameters(train, test)

      #Check for constraints
      stopText <- ""

      if((hasMissingValues(train)) || (hasMissingValues(test))){
        stopText <- paste0(stopText, "Dataset has missing values and the algorithm does not accept it.\n")
      }

      if(stopText != ""){
        stop(stopText)
      }

      self$C <- C
      self$toleranceParameter <- toleranceParameter
      self$epsilon <- epsilon
      self$RBFKernel_gamma <- RBFKernel_gamma
      self$normalized_PolyKernel_exponent <- normalized_PolyKernel_exponent
      self$normalized_PolyKernel_useLowerOrder <- normalized_PolyKernel_useLowerOrder
      self$PukKernel_omega <- PukKernel_omega
      self$PukKernel_sigma <- PukKernel_sigma
      self$StringKernel_lambda <- StringKernel_lambda
      self$StringKernel_subsequenceLength <- StringKernel_subsequenceLength
      self$StringKernel_maxSubsequenceLength <- StringKernel_maxSubsequenceLength
      self$StringKernel_normalize <- StringKernel_normalize

      if((tolower(StringKernel_pruning) == "none") || (tolower(StringKernel_pruning) == "lambda")){
        self$StringKernel_pruning <- StringKernel_pruning
      }
      else{
        #Default value
        self$StringKernel_pruning <- "None"
      }

      if((tolower(KernelType) == "rbfkernel") || (tolower(KernelType) == "polykernel") || (tolower(KernelType) == "normalizedpolykernel") || (tolower(KernelType) == "puk") || (tolower(KernelType) == "stringkernel")){
        self$KernelType <- KernelType
      }
      else{
        #Default value
        self$StringKernel_pruning <- "PolyKernel"
      }

      self$FitLogisticModels <- FitLogisticModels
      self$ConvertNominalAttributesToBinary <- ConvertNominalAttributesToBinary

      if((tolower(PreprocessType) == "normalize") || (tolower(PreprocessType) == "standarize") || (tolower(PreprocessType) == "none")){
        self$PreprocessType <- PreprocessType
      }
      else{
        #Default value
        self$PreprocessType <- "Normalize"
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
    jarName = "SMO.jar",

    #algorithm name
    algorithmName = "SMO-C",

    #String with algorithm name
    algorithmString = "SMO",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "C = ", self$C, "\n")
      text <- paste0(text, "toleranceParameter = ", self$toleranceParameter, "\n")
      text <- paste0(text, "epsilon = ", self$epsilon, "\n")
      text <- paste0(text, "RBFKernel_gamma = ", self$RBFKernel_gamma, "\n")
      text <- paste0(text, "-Normalized-PolyKernel_exponent = ", self$normalized_PolyKernel_exponent, "\n")
      text <- paste0(text, "-Normalized-PolyKernel_useLowerOrder = ", self$normalized_PolyKernel_useLowerOrder, "\n")
      text <- paste0(text, "PukKernel_omega = ", self$PukKernel_omega, "\n")
      text <- paste0(text, "PukKernel_sigma = ", self$PukKernel_sigma, "\n")
      text <- paste0(text, "StringKernel_lambda = ", self$StringKernel_lambda, "\n")
      text <- paste0(text, "StringKernel_subsequenceLength = ", self$StringKernel_subsequenceLength, "\n")
      text <- paste0(text, "StringKernel_maxSubsequenceLength = ", self$StringKernel_maxSubsequenceLength, "\n")
      text <- paste0(text, "StringKernel_normalize = ", self$StringKernel_normalize, "\n")
      text <- paste0(text, "StringKernel_pruning = ", self$StringKernel_pruning, "\n")
      text <- paste0(text, "KERNELtype = ", self$KernelType, "\n")
      text <- paste0(text, "FitLogisticModels = ", self$FitLogisticModels, "\n")
      text <- paste0(text, "ConvertNominalAttributesToBinary = ", self$ConvertNominalAttributesToBinary, "\n")
      text <- paste0(text, "PreprocessType = ", self$PreprocessType, "\n")

      return(text)

    }

  )
)

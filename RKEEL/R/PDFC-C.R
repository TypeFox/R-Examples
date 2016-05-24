#Class implementing a Classification Algorithm
  #Implements the PDFC-C KEEL classification algorithm

PDFC_C <- function(train, test, C=100.0, d=0.25, tolerance=0.001, epsilon=1.0E-12, PDRFtype="Gaussian", nominal_to_binary=TRUE, preprocess_type="Normalize", seed=-1){
  alg <- RKEEL::R6_PDFC_C$new()
  alg$setParameters(train, test, C, d, tolerance, epsilon, PDRFtype, nominal_to_binary, preprocess_type, seed)
  return (alg)
}

R6_PDFC_C <- R6::R6Class("R6_PDFC_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #C
    C = 100.0,

    #d
    d = 0.25,

    #Tolerance parameter
    tolerance = 0.001,

    #epsilon
    epsilon = 1.0E-12,

    #PDRF type
    PDRFtype = "Gaussian",

    #Convert nominal attributes to binary?
    nominal_to_binary = TRUE,

    #Preprocess type
    preprocess_type = "Normalize",

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, C=100.0, d=0.25, tolerance=0.001,
                          epsilon=1.0E-12, PDRFtype="Gaussian",
                          nominal_to_binary=TRUE, preprocess_type="Normalize",
                          seed=-1){

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
      self$d <- d
      self$tolerance <- tolerance
      self$epsilon <- epsilon


      if((tolower(PDRFtype) == "symmetrictriangle")){
        self$PDRFtype <- "SymmetricTriangle"
      }
      else if((tolower(PDRFtype) == "gaussian")){
        self$PDRFtype <- "Gaussian"
      }
      else if((tolower(PDRFtype) == "cauchy")){
        self$PDRFtype <- "Cauchy"
      }
      else if((tolower(PDRFtype) == "laplace")){
        self$PDRFtype <- "Laplace"
      }
      else if((tolower(PDRFtype) == "hyperbolicsecant")){
        self$PDRFtype <- "HyperbolicSecant"
      }
      else if((tolower(PDRFtype) == "squaredsinc")){
        self$PDRFtype <- "SquaredSinc"
      }
      else{
        #Default value
        self$PDRFtype <- "Gaussian"
      }

      self$nominal_to_binary <- nominal_to_binary

      if((tolower(preprocess_type) == "normalize")){
        self$preprocess_type <- "Normalize"
      }
      else if((tolower(preprocess_type) == "standardize")){
        self$preprocess_type <- "Standardize"
      }
      else if((tolower(preprocess_type) == "none")){
        self$preprocess_type <- "None"
      }
      else{
        #Default value
        self$preprocess_type <- "Normalize"
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
    jarName = "PDFC.jar",

    #algorithm name
    algorithmName = "PDFC-C",

    #String with algorithm name
    algorithmString = "Positive Definite Fuzzy Classifier",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "C = ", self$C, "\n")
      text <- paste0(text, "d = ", self$d, "\n")
      text <- paste0(text, "toleranceParameter = ", self$tolerance, "\n")
      text <- paste0(text, "epsilon = ", self$epsilon, "\n")
      text <- paste0(text, "PDRFtype = ", self$PDRFtype, "\n")
      if(self$nominal_to_binary){
        text <- paste0(text, "ConvertNominalAttributesToBinary = True", "\n")
      }
      else{
        text <- paste0(text, "ConvertNominalAttributesToBinary = False", "\n")
      }
      text <- paste0(text, "PreprocessType = ", self$preprocess_type, "\n")

      return(text)

    }

  )
)

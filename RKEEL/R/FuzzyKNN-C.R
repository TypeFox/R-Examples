#Class implementing a Classification Algorithm
  #Implements the FuzzyKNN-C KEEL classification algorithm

FuzzyKNN_C <- function(train, test, k=3, M=2.0, initialization="CRISP", init_k=3){
  alg <- RKEEL::R6_FuzzyKNN_C$new()
  alg$setParameters(train, test, k, M, initialization, init_k)
  return (alg)
}

R6_FuzzyKNN_C <- R6::R6Class("R6_FuzzyKNN_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #k
    k = 3,

    #M
    M = 2.0,

    #initialization
    initialization = "CRISP",

    #k value (init)
    init_k = 3,

    #Public functions

    #Initialize function
    setParameters = function(train, test, k=3, M=2.0, initialization="CRISP",
                          init_k=3){

      super$setParameters(train, test)

      self$k <- k
      self$M <- M

      if((tolower(initialization) == "crisp") || (tolower(initialization) == "class-mean") || (tolower(initialization) == "knn")){
        self$initialization <- toupper(initialization)
      }
      else{
        #Default value of distance type
        self$distance <- "CRISP"
      }

      self$init_k <- init_k

    }

  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "FuzzyKNN.jar",

    #algorithm name
    algorithmName = "FuzzyKNN-C",

    #String with algorithm name
    algorithmString = "FuzzyKNN",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "K Value = ", self$k, "\n")
      text <- paste0(text, "M = ", self$M, "\n")
      text <- paste0(text, "Initialization = ", self$initialization, "\n")
      text <- paste0(text, "K Value (init) = ", self$init_k, "\n")

      return(text)

    }

  )
)

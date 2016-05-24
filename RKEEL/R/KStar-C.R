#Class implementing a Classification Algorithm
  #Implements the KStar-C KEEL classification algorithm

KStar_C <- function(train, test, selection_method="Fixed", blend=0.2, seed=-1){
  alg <- RKEEL::R6_KStar_C$new()
  alg$setParameters(train, test, selection_method, blend, seed)
  return (alg)
}

R6_KStar_C <- R6::R6Class("R6_KStar_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #Selection method
    selection_method = "Fixed",

    #Blend
    blend = 0.2,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, selection_method="Fixed", blend=0.2,
                          seed=-1){

      super$setParameters(train, test)

      if(tolower(selection_method) == "fixed"){
        self$selection_method = "Fixed"
      }
      else if(tolower(selection_method) == "random"){
        self$selection_method = "Random"
      }
      else{
        #Default value
        self$selection_method = "Fixed"
      }

      self$blend <- blend

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
    jarName = "KStar.jar",

    #algorithm name
    algorithmName = "KStar-C",

    #String with algorithm name
    algorithmString = "K Star Classifier",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "Selection method = ", self$selection_method, "\n")
      text <- paste0(text, "Blend = ", self$blend, "\n")

      return(text)

    }
  )
)

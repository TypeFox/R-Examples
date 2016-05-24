#Class implementing a Classification Algorithm
  #Implements the BSE-C KEEL classification algorithm

BSE_C <- function(train, test, k=1, distance="Euclidean"){
  alg <- RKEEL::R6_BSE_C$new()
  alg$setParameters(train, test, k, distance)
  return (alg)
}

R6_BSE_C <- R6::R6Class("R6_BSE_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #k: num of nearest neighbours
    k = 1,

    #distance type
    distance = "Euclidean",


    #Public functions

    #Initialize function
    setParameters = function(train, test, k=1, distance="Euclidean"){

      super$setParameters(train, test)

      self$k <- k

      if((tolower(distance) == "euclidean") || (tolower(distance) == "hvdm")){
        self$distance <- distance
      }
      else{
        #Default value of distance type
        self$distance <- "Euclidean"
      }

    }

  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "IS-BSE.jar",

    #algorithm name
    algorithmName = "BSE-C",

    #String with algorithm name
    algorithmString = "Backward Sequential Edition",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "K Value = ", self$k, "\n")
      text <- paste0(text, "Distance Function = ", self$distance, "\n")

      return(text)

    }

  )
)

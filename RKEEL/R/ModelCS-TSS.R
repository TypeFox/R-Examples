#Class implementing a Preprocess Algorithm
  #Implements the ModelCS-TSS KEEL preprocess algorithm

ModelCS_TSS <- function(train, test, k=3, distance="Euclidean"){
  alg <- RKEEL::R6_ModelCS_TSS$new()
  alg$setParameters(train, test, k, distance)
  return (alg)
}

R6_ModelCS_TSS <- R6::R6Class("R6_ModelCS_TSS",

  inherit = PreprocessAlgorithm,

  public = list(

    #Public properties

    #Number of Neighbors
    k = 3,

    #distance function
    distance = "Euclidean",


    #Public functions

    #Initialize function
    setParameters = function(train, test, k=3, distance="Euclidean"){

      super$setParameters(train, test)

      self$k <- k

      if((tolower(distance) == "euclidean") || (tolower(distance) == "hvdm")){
        self$distance <- distance
      }
      else{
        #Default value of distance function
        self$distance <- "Euclidean"
      }

    }


  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "IS-ModelCS.jar",

    #algorithm name
    algorithmName = "ModelCS-TSS",

    #String with algorithm name
    algorithmString = "Model Class Selection",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "Number of Neighbors = ", self$k, "\n")
      text <- paste0(text, "Distance Function = ", self$distance, "\n")

      return(text)

    }
  )
)

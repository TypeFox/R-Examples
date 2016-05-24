#Class implementing a Classification Algorithm
  #Implements the CNN-C KEEL classification algorithm

CNN_C <- function(train, test, k=1, distance="Euclidean", seed=-1){
  alg <- RKEEL::R6_CNN_C$new()
  alg$setParameters(train, test, k, distance, seed)
  return (alg)
}

R6_CNN_C <- R6::R6Class("R6_CNN_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #k: num of nearest neighbours
    k = 1,

    #distance type
    distance = "Euclidean",

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, k=1, distance="Euclidean", seed=-1){

      super$setParameters(train, test)

      self$k <- k

      if((tolower(distance) == "euclidean") || (tolower(distance) == "hvdm")){
        self$distance <- distance
      }
      else{
        #Default value of distance type
        self$distance <- "Euclidean"
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
    jarName = "IS-CNN.jar",

    #algorithm name
    algorithmName = "CNN-C",

    #String with algorithm name
    algorithmString = "Condensed Nearest Neighbor",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "K Value = ", self$k, "\n")
      text <- paste0(text, "Distance Function = ", self$distance, "\n")

      return(text)

    }

  )
)

#Class implementing a Classification Algorithm
  #Implements the C45_Binarization-C KEEL classification algorithm

C45Binarization_C <- function(train, test, pruned=TRUE, confidence=0.25, instancesPerLeaf=2, binarization="OVO", scoreFunction="WEIGHTED", bts=0.05){
  alg <- RKEEL::R6_C45Binarization_C$new()
  alg$setParameters(train, test, pruned, confidence, instancesPerLeaf, binarization, scoreFunction, bts)
  return (alg)
}

R6_C45Binarization_C <- R6::R6Class("R6_C45Binarization_C",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public properties

    #pruned
    pruned = TRUE,

    #confidence
    confidence = 0.25,

    #instances per leaf
    instancesPerLeaf = 2,

    #Binarization
    binarization = "OVO",

    #Score function
    scoreFunction = "WEIGHTED",

    #BTS Param (if applied)
    bts = 0.05,

    #Public functions

    #Initialize function
    setParameters = function(train, test, pruned=TRUE, confidence=0.25,
                          instancesPerLeaf=2, binarization="OVO",
                          scoreFunction="WEIGHTED", bts=0.05){

      super$setParameters(train, test)

      self$pruned <- pruned

      self$confidence <- confidence

      self$instancesPerLeaf <- 2

      if((tolower(binarization) == "ovo") || (tolower(binarization) == "ova")){
        self$binarization <- toupper(binarization)
      }
      else{
        #Default value of training method
        self$trainMethod <- "OVO"
      }

      if((tolower(scoreFunction) == "vote") || (tolower(scoreFunction) == "weighted") ||
         (tolower(scoreFunction) == "ddag") || (tolower(scoreFunction) == "pc") ||
         (tolower(scoreFunction) == "nd") || (tolower(scoreFunction) == "lvpc") ||
         (tolower(scoreFunction) == "bts") || (tolower(scoreFunction) == "wujmlr") ||
         (tolower(scoreFunction) == "nesting") || (tolower(scoreFunction) == "dynovo")){

        if(tolower(scoreFunction) == "wujmlr"){
          self$scoreFunction <- "WuJMLR"
        }
        else if(tolower(scoreFunction) == "dynovo") {
          self$scoreFunction <- "DynOVO"
        }
        else {
          self$scoreFunction <- toupper(scoreFunction)
        }
      }
      else{
        #Default value of training method
        self$scoreFunction <- "WEIGHTED"
      }

      self$bts <- bts

    }

  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "C45_Binarization-C.jar",

    #algorithm name
    algorithmName = "C45_Binarization-C",

    #String with algorithm name
    algorithmString = "C4.5-OVO",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "pruned = ", self$pruned, "\n")
      text <- paste0(text, "confidence = ", self$confidence, "\n")
      text <- paste0(text, "isntancesPerLeaf = ", self$instancesPerLeaf, "\n")
      text <- paste0(text, "Binarization = ", self$binarization, "\n")
      text <- paste0(text, "Score Function = ", self$scoreFunction, "\n")
      text <- paste0(text, "BTS_Param (if applied) = ", self$bts, "\n")

      return(text)

    }
  )
)

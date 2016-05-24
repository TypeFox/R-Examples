#Class implementing a Preprocess Algorithm
  #Implements the LVF-IEP-FS KEEL preprocess algorithm

LVF_IEP_FS <- function(train, test, paramKNN=1, maxLoops=770, inconAllow=0, seed=-1){
  alg <- RKEEL::R6_LVF_IEP_FS$new()
  alg$setParameters(train, test, paramKNN, maxLoops, inconAllow, seed)
  return (alg)
}

R6_LVF_IEP_FS <- R6::R6Class("R6_LVF_IEP_FS",

  inherit = PreprocessAlgorithm,

  public = list(

    #Public properties

    #paramKNN
    paramKNN = 1,

    #Max loops
    maxLoops = 770,

    #inconAllow
    inconAllow = 0,

    #seed
    seed = -1,

    #Public functions

    #Initialize function
    setParameters = function(train, test, paramKNN=1, maxLoops=770, inconAllow=0,
                          seed=-1){

      super$setParameters(train, test)

      stopText <- ""

      if(is.null(test)){
        if(hasMissingValues(train)){
          stopText <- paste0(stopText, "Dataset has missing values and the algorithm does not accept it.\n")
        }
        if(hasContinuousData(train)){
          stopText <- paste0(stopText, "Dataset has continuous data and the algorithm does not accept it.\n")
        }
      }
      else{
        if((hasMissingValues(train)) || (hasMissingValues(test))){
          stopText <- paste0(stopText, "Dataset has missing values and the algorithm does not accept it.\n")
        }
        if((hasContinuousData(train)) || (hasContinuousData(test))){
          stopText <- paste0(stopText, "Dataset has continuous data and the algorithm does not accept it.\n")
        }
      }

      if(stopText != ""){
        stop(stopText)
      }

      self$paramKNN <- paramKNN
      self$maxLoops <- maxLoops
      self$inconAllow <- inconAllow

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
    jarName = "FS-LVF-IEP.jar",

    #algorithm name
    algorithmName = "LVF-IEP-FS",

    #String with algorithm name
    algorithmString = "Las Vegas Filter (IEP)",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")
      text <- paste0(text, "paramKNN = ", self$paramKNN, "\n")
      text <- paste0(text, "maxLoops = ", self$maxLoops, "\n")
      text <- paste0(text, "inconAllow = ", self$inconAllow, "\n")

      return(text)

    }
  )
)

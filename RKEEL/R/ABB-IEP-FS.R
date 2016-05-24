#Class implementing a Preprocess Algorithm
  #Implements the ABB-IEP-FS KEEL preprocess algorithm

ABB_IEP_FS <- function(train, test, seed=-1){
  alg <- RKEEL::R6_ABB_IEP_FS$new()
  alg$setParameters(train, test, seed)
  return (alg)
}

R6_ABB_IEP_FS <- R6::R6Class("R6_ABB_IEP_FS",

  inherit = PreprocessAlgorithm,

  public = list(

    #Public properties

    #seed
    seed = -1,

    #Public functions

    #Initialize function
    setParameters = function(train, test, seed=-1){

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
    jarName = "FS-ABB-IEP.jar",

    #algorithm name
    algorithmName = "ABB-IEP-FS",

    #String with algorithm name
    algorithmString = "Automatic Branch and Bound (IEP)",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "seed = ", self$seed, "\n")

      return(text)

    }
  )
)

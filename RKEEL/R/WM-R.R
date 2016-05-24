#Class implementing a Regression Algorithm
  #Implements the WM-R KEEL regression algorithm

WM_R <- function(train, test, numlabels=5, KB=FALSE){
  alg <- RKEEL::R6_WM_R$new()
  alg$setParameters(train, test, numlabels, KB)
  return (alg)
}

R6_WM_R <- R6::R6Class("R6_WM_R",

  inherit = RegressionAlgorithm,

  public = list(

    #Public properties

    #number of labels
    numlabels = 5,

    #KB Output File Format with Weight values to 1? (0/1)
    KB = FALSE,

    #seed
    seed = -1,


    #Public functions

    #Initialize function
    setParameters = function(train, test, numlabels=5, KB=FALSE){

      super$setParameters(train, test)

      self$numlabels <- numlabels
      self$KB <- KB
    }

  ),

  private = list(

    #Private properties

    #jar Filename
    jarName = "Mam-WM.jar",

    #algorithm name
    algorithmName = "WM-R",

    #String with algorithm name
    algorithmString = "Fuzzy Rule Learning, Wang-Mendel Algorithm",


    #Private functions

    #Get the text with the parameters for the config file
    getParametersText = function(){

      text <- ""
      text <- paste0(text, "Number of labels = ", self$numlabels, "\n")
      if(self$KB == FALSE){
        text <- paste0(text, "KB Output File Format with Weight values to 1? (0/1) = 0", "\n")
      }
      else{
        text <- paste0(text, "KB Output File Format with Weight values to 1? (0/1) = 1", "\n")
      }

      return(text)

    },

    #Create config files function
    writeKeelConfig = function(){

      #Ficheros de entrada
      inputDataString = paste0("\"../datasets/", private$dataName, "/", private$trainFilename, "\" \"../datasets/", private$dataName, "/", private$trainFilename, "\" \"../datasets/", private$dataName, "/", private$testFilename, "\"")
      #Ficheros de salida
      #dataName <- strsplit(trainData, "/")[[1]][1]
      outputDataString = paste0("\"../results/", private$algorithmName, ".",private$dataName, "/result0", ".tra\" \"../results/", private$algorithmName, ".", private$dataName, "/result0", ".tst\" \"../results/", private$algorithmName, ".", private$dataName, "/result0", "e0.txt\" \"../results/", private$algorithmName, ".", private$dataName, "/result0", "e1.txt\"")

      text <- ""
      text <- paste0(text, "algorithm = ", private$algorithmString)
      text <- paste0(text, "\ninputData = ", inputDataString)
      text <- paste0(text, "\noutputData = ", outputDataString)
      text <- paste0(text, "\n\n")
      text <- paste0(text, private$getParametersText())
      text <- paste0(text, "\n")

      writeLines(text, con = paste0(private$mainPath, "/scripts/", private$algorithmName, "/", private$dataName, "/", "config0s0", ".txt"))

    }

  )
)

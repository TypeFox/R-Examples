#Class implementing a Classification Algorithm
#Implements specific functions of an Associative classification algorithm

AssociativeClassificationAlgorithm <- R6::R6Class("AssociativeClassificationAlgorithm",

  inherit = ClassificationAlgorithm,

  public = list(

    #Public functions

    #Initialize function
    setParameters = function(train, test){
      super$setParameters(train, test)
    },

    run = function(){
      super$run()
    },

    #Read algorithm output files
    readOutputs = function(trafile, tstfile, txtfile){

      #Check prediction files
      if((! file.exists(trafile)) || (! file.exists(tstfile))){
        cat(trafile, "\n", sep="")
        cat(tstfile, "\n", sep="")
        stop("Results were not stored. Execution error.")
      }

      #Read text output file only if it exists
      if(! file.exists(txtfile)) {
        self$output <- NULL
      }
      else {
        #Read output from txt file
        self$output <- readLines(txtfile, warn = FALSE)
      }

      txtfile2 <- gsub("result0e0", "result0e1", txtfile)
      #Read extra text output file only if it exists
      if(file.exists(txtfile2)) {
        #Read output from txt file
        self$output <- paste0(self$output, "\n\n", "*----------*", "\n\n", readLines(txtfile2, warn = FALSE))
      }

      #Read predictions
      self$trainPredictions <- private$readPredictionOutputs(trafile)
      self$testPredictions <- private$readPredictionOutputs(tstfile)
    }

  ),

  private = list(

    #Private functions

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

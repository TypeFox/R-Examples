#Class implementing the KeelAlgorithm
  #Implements the common functions of a KEEL preprocess algorithm

PreprocessAlgorithm <- R6::R6Class("PreprocessAlgorithm",

  inherit = KeelAlgorithm,

  public = list(

    #Public properties

    #Algorithm outputs
    preprocessed_train = NULL,
    preprocessed_test = NULL,
    output = NULL,

    #Public functions

    #Initialize function
    setParameters = function(train, test=NULL){

      private$trainFilename <- "train.dat"
      private$trainDataset <- train

      if(is.null(test)) {
        private$testFilename <- "train.dat"
        private$testDataset <- NULL
      }
      else {
        private$testFilename <- "test.dat"
        private$testDataset <- test
      }

      #super$initialize()

      #Test jar file
      #if(! file.exists(paste0(private$exePath, private$jarName))){
      if(! file.exists(system.file("exe", private$jarName, package = "RKEELjars"))){
        stop(paste0(private$jarName, " doesn't exist under the defined path. Installation error."))
      }

      private$dataName <- "data"
    },

    run = function(){

      #Use tryCatch() to remove experiment folders even it there are errors
      tryCatch({
        #Experiment folder
        #expPath <- gsub("//", "/", system.file("exp", "", package="RKEEL"))

        expPath <- gsub("\\\\", "/", tempdir())

        if(substr(expPath, nchar(expPath), nchar(expPath)) != "/"){
          expPath <- paste0(expPath, "/")
        }

        private$mainPath <- paste0(expPath, "experiment_", gsub(" ", "_", gsub(":", "-", toString(Sys.time()))), sample(1:10000, 1))
        private$generateExperimentDir(private$mainPath)

        #Copy dataset folder
        #Create dataset folder
        dir.create(paste0(private$mainPath, "/datasets/", private$dataName))
        dir.create(paste0(private$mainPath, "/datasets/", private$algorithmName, ".", private$dataName))
        #Write dataset files
        if(is.null(private$testDataset)){
          writeDatFromDataframe(private$trainDataset, paste0(private$mainPath, "/datasets/", private$dataName, "/", private$trainFilename))
        }
        else {
          writeDatFromDataframes(private$trainDataset, private$testDataset, paste0(private$mainPath, "/datasets/", private$dataName, "/", private$trainFilename), paste0(private$mainPath, "/datasets/", private$dataName, "/", private$testFilename))
        }

        #Copy algorithm exe
        file.copy(system.file("exe", private$jarName, package = "RKEELjars"), paste0(private$mainPath, "/exe/", private$jarName))

        #Create results dir
        #dir.create(paste0(private$mainPath, "/results/", private$algorithmName, ".", private$dataName))

        #Create .xml experiment file
        private$writeKeelXML()

        #Create config files
        dir.create(paste0(private$mainPath, "/scripts/", private$algorithmName))
        dir.create(paste0(private$mainPath, "/scripts/", private$algorithmName, "/", private$dataName))
        private$writeKeelConfig()

        #Change work directory to execute .jar
        wdPath <- getwd()
        setwd(paste0(private$mainPath, "/scripts/"))
        if(grepl("windows", tolower(Sys.info()["sysname"]))) {
          system("java -jar RunKeel.jar", show.output.on.console = FALSE)
        }
        else {
          system("java -jar RunKeel.jar", ignore.stdout = TRUE)
        }
        setwd(wdPath)

        #read outputs
        of1 <- paste0(private$mainPath, "/datasets/", private$algorithmName, ".",private$dataName, "/train.dat")
        if(is.null(private$testDataset)){
          of2 <- NULL
        }
        else {
          of2 <- paste0(private$mainPath, "/datasets/", private$algorithmName, ".",private$dataName, "/test.dat")
        }
        of3 <- paste0(private$mainPath, "/datasets/", private$algorithmName, ".",private$dataName, "/results.txt")
        private$readOutputs(of1, of2, of3)


        cat(paste0("Algorithm ",  class(self)[1], " executed successfully", sep="\n"))

      }, error = function(err) {
        #Error
        cat(paste0("Error! ",err))
      }, finally = {
        #Remove data files and Keel experiment folder
        unlink(paste0(private$dataPath, private$trainFilename))
        unlink(paste0(private$dataPath, private$testFilename))
        unlink(private$mainPath, recursive = TRUE)
      })

    },


    #Print object
    print = function() {

      if(class(self)[1] == "PreprocessAlgorithm") {
        cat("Keel Preprocess Algorithm Object", sep="\n")
        cat("Need to have an object of an algorithm implementing this class")
      }
      else {
        cat("-----------------------------------", "\n", sep="")
        cat(private$algorithmString, "\n", sep="")
        cat("-----------------------------------", "\n", sep="")
        cat(private$getParametersText(), sep="")
        cat("-----------------------------------", "\n", sep="")
      }

    }

  ),

  private = list(

    #train filename
    trainFilename = NULL,

    #test filename
    testFilename = NULL,

    #train dataset
    trainDataset = NULL,

    #test dataset
    testDataset = NULL,

    #dataset name
    dataName = NULL,

    #experiment main path
    mainPath = NULL,

    #Create XML function
    writeKeelXML = function(){

      nodeExecution <- XML::xmlNode("execution")

      nodeSentence <- XML::xmlNode("sentence")
      nodeSentence <- XML::addChildren(nodeSentence, kids = list(XML::xmlNode("command", "java")))
      nodeSentence <- XML::addChildren(nodeSentence, kids = list(XML::xmlNode("option", "-Xmx512000000")))
      nodeSentence <- XML::addChildren(nodeSentence, kids = list(XML::xmlNode("option", "-jar")))
      nodeSentence <- XML::addChildren(nodeSentence, kids = list(XML::xmlNode("executableFile", paste0("../exe/", private$jarName))))
      nodeSentence <- XML::addChildren(nodeSentence, kids = list(XML::xmlNode("scriptFile", paste0("./", private$algorithmName, "/", private$dataName, "/config0s0", ".txt"))))
      nodeExecution <- XML::addChildren(nodeExecution, kids = list(nodeSentence))

      XML::saveXML(nodeExecution, paste0(private$mainPath, "/scripts/", "RunKeel.xml"))

    },


    #Create config files function
    writeKeelConfig = function(){

      #Ficheros de entrada
      if(is.null(private$testDataset)) {
        inputDataString = paste0("\"../datasets/", private$dataName, "/", private$trainFilename, "\" \"../datasets/", private$dataName, "/", private$trainFilename, "\"")
      }
      else {
        inputDataString = paste0("\"../datasets/", private$dataName, "/", private$trainFilename, "\" \"../datasets/", private$dataName, "/", private$testFilename, "\"")
      }
      #Ficheros de salida
      #dataName <- strsplit(trainData, "/")[[1]][1]
      outputDataString = paste0("\"../datasets/", private$algorithmName, ".",private$dataName, "/train.dat\" \"../datasets/", private$algorithmName, ".", private$dataName, "/test.dat\" ", "\"../datasets/", private$algorithmName, ".",private$dataName, "/results.txt\"")

      text <- ""
      text <- paste0(text, "algorithm = ", private$algorithmString)
      text <- paste0(text, "\ninputData = ", inputDataString)
      text <- paste0(text, "\noutputData = ", outputDataString)
      text <- paste0(text, "\n\n")
      text <- paste0(text, private$getParametersText())
      text <- paste0(text, "\n")

      writeLines(text, con = paste0(private$mainPath, "/scripts/", private$algorithmName, "/", private$dataName, "/", "config0s0", ".txt"))

    },

    #Read algorithm output files
    readOutputs = function(trafile, tstfile, txtfile){

      #Check prediction files
      if(is.null(private$testDataset)){
        if(! file.exists(trafile)){
          cat(getwd(), "\n", sep="")
          cat(trafile, "\n", sep="")
          stop("Results were not stored. Execution error.")
        }
      }
      else{
        if((! file.exists(trafile)) || (! file.exists(tstfile))){
          cat(getwd(), "\n", sep="")
          cat(trafile, "\n", sep="")
          cat(tstfile, "\n", sep="")
          stop("Results were not stored. Execution error.")
        }
      }

      self$preprocessed_train <- read.keel(trafile)
      if(is.null(private$testDataset)){
        self$preprocessed_test <- NULL
      }
      else{
        self$preprocessed_test <- read.keel(tstfile)
      }

      #Read text output file only if it exists
      if(! file.exists(txtfile)) {
        self$output <- NULL
      }
      else {
        #Read file
        self$output <- readLines(txtfile, warn = FALSE)
      }
    }

  )
)

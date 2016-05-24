#Class that defines a KEEL Algorithm
  #Implements the common functions of a KEEL algorithm

KeelAlgorithm <- R6::R6Class("KeelAlgorithm",

  public = list(

    #Public properties


    #Public functions

    #Initialize
    initialize = function(){
      jarPath <- RKEELjars::getJarPath()

      if(substr(private$jarPath, nchar(private$jarPath), nchar(private$jarPath)) != "/"){
        private$jarPath <- paste0(private$jarPath, "/")
      }
      if(substr(private$exePath, nchar(private$exePath), nchar(private$exePath)) != "/"){
        private$exePath <- paste0(private$exePath, "/")
      }
      if(substr(private$dataPath, nchar(private$dataPath), nchar(private$dataPath)) != "/"){
        private$dataPath <- paste0(private$dataPath, "/")
      }

      #Test paths
      #cat(paste0("HOLA:", system.file("exe", "RunKeel.jar", package = "RKEELjars")))

      if(! file.exists(system.file("exe", "RunKeel.jar", package = "RKEELjars"))){
        stop("RunKeel.jar doesn't exist under the defined path. Installation error.")
      }

      #if(! file.exists(paste0(private$jarPath, "RunKeel.jar"))) {
      # cat(paste0("jarPath: ", private$jarPath), sep="\n")
      # cat(paste0(private$jarPath, "RunKeel.jar"))
      # stop("RunKeel.jar doesn't exist under the defined path. Installation error.")
      #}
    },

    #Execute algorithm
    run = function(){
      #Implement in each algorithm type class
    },

    #Print object
    print = function(...) {
      cat("Keel Algorithm Object", sep="\n")
      cat("Need to have an object of an algorithm implementing this class")
    }

  ),

  private = list(

    #Private properties

    #RunKeel.jar path
    jarPath = RKEELjars::getJarPath(),

    #exes path
    exePath = RKEELjars::getExePath(),

    #dataset path
    dataPath = RKEELdata::getDataPath(),


    #Private Functions

    #Create XML function
    writeKeelXML = function(mainPath){
      #Implement in each algorithm class
    },

    #Create config files function
    writeKeelConfig = function(mainPath){
      #Implement in each algorithm class
    },

    #Generate experiment directory
    generateExperimentDir = function(mainPath){
      dir.create(mainPath)
      dir.create(paste0(mainPath, "/datasets"))
      dir.create(paste0(mainPath, "/exe"))
      dir.create(paste0(mainPath, "/results"))
      dir.create(paste0(mainPath, "/scripts"))
      #file.copy(paste0(private$jarPath, "RunKeel.jar"), paste0(mainPath, "/scripts/RunKeel.jar"))
      file.copy(system.file("exe", "RunKeel.jar", package = "RKEELjars"), paste0(mainPath, "/scripts/RunKeel.jar"))
    }

  )
)

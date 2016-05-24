##================================================================================
## This file is part of the rrepast package - R/Repast interface API
## 
## (C)2015 Antonio Prestes Garcia <@>
## For license terms see DESCRIPTION and/or LICENSE
##
## $Id$
##================================================================================


# ------------------------------------------------------------
# .onLoad, Hook for loading package namespace
# 
# ------------------------------------------------------------
.onLoad<- function(libname, pkgname) {
  #packageStartupMessage("R/Repast: Integrating Repast Models into R\n")
  assign("pkg.globals", new.env(), envir=parent.env(environment()))
  
  # Internal variables
  assign("pkg.basedir", NA, pkg.globals)
  assign("pkg.modeldir", NA, pkg.globals)
  assign("pkg.scenariodir", NA, pkg.globals)
  assign("pkg.modellibdir", NA, pkg.globals)
  assign("pkg.id", NA, pkg.globals)
  
  # global simulation results
  assign("pkg.parameters", data.frame(), pkg.globals)
  assign("pkg.results", data.frame(), pkg.globals)
  
  # progress bar internals
  assign("pkg.progressbar", NULL, pkg.globals)
  assign("pkg.progressbar.enabled", FALSE, pkg.globals)

  # default values for model
  assign("pkg.outputdir",paste0(Sys.getenv("TMP"),"/rrepast-deployment/"), pkg.globals)
  assign("pkg.repastlibdir", "/repast.simphony/", pkg.globals)
  assign("pkg.java.parameters","-server -Xms512m -Xmx1024m", pkg.globals)
  
  # default key for Repast random seed 
  assign("pkg.randomSeed","randomSeed", pkg.globals)
  
  # The Random Seed. You may want to change this.
  set.seed(exp(1)*10^6)
  
  # Define funcions which are not present in ond R versions
  compatibility()
}



# ----- internal functions 

# Define some required functions when not available from current R version
compatibility<- function() {
  if(getRversion() <= "3.1") {
    # dir.exists function is only available from R 3.2
    f<- function(d) {
      cat("my dir.exists()")
      v<- file.info(d)$isdir
      return(ifelse(is.na(v), FALSE, v))
    }
    
    # This is a trick to get the reference to the current environment
    e<- as.environment(environment(enginejar))
    assign("dir.exists", f, e, immediate = FALSE)
  }
}

# Returns the wrapper classes jar file location. 
enginejar<- function() {

  # Try to guess the rrepast-engine.jar location  
  for(p in .libPaths()) {
    f0<- paste0(p,"/rrepast/java/rrepast-engine.jar")
    f1<- paste0(p,"/rrepast/inst/java/rrepast-engine.jar")
    if(file.exists(f0)) {
      f<- f0
      break
    } else if(file.exists(f1)) {
        f<- f1
        break
    }
  }
  return(f)
}

# Return the name of repast engine class name. 
engineclazz<- function() {
  return("org.haldane.rrepast.RepastEngine")
}

# Returns the repast simphony library dir
simphonylib<- function() {
  b<- get("pkg.basedir", pkg.globals)
  l<- get("pkg.repastlibdir", pkg.globals)
  return(paste0(b,l))
}

# Configure all model directories based on default installation values
configModelDirs<- function(s) {
  d<- basename(s)
  setBaseDir(s)
  setModelDir(paste0(paste0(s,"/"),d))
  setScenarioDir(paste0(paste0(paste0(getModelDir(),"/"),d),".rs"))
  setModelLibDir(paste0(getModelDir(),"/lib"))
}

# Setters and Getters ----------

#' @title Sets the model name
#' @description  Set the name of the model currently instantiated.
#' 
#' @param s The model name
#' 
#' @export
setId<- function(s) {
  assign("pkg.id", s, pkg.globals)      
}

#' @title Gets the model name
#' @description  Provides the name of the model currently instantiated.
#' 
#' @export
getId<- function() {
  return(get("pkg.id", pkg.globals))
}

#' @title Sets Repast randomSeed name
#' @description Configures a non-default value for Repast randomSeed
#' parameter name.
#' 
#' @param k The string with an alternative name for randomSeed 
#' 
#' @export
setKeyRandom<- function(k){
  assign("pkg.randomSeed",k, pkg.globals)
}

#' @title Gets Repast randomSeed name
#' @description Returns the Repast randomSeed parameter name.
#' 
#' @return A string value holding the randomSeed name.
#' 
#' @export
getKeyRandom<- function() {
  return(get("pkg.randomSeed", pkg.globals))
}

#' @title Sets output directory
#' 
#' @description Configure the desired directoy to save model 
#' output data.
#' 
#' @param s The full path for output directory
#' 
#' @export
setOutputDir<- function(s) {
  assign("pkg.outputdir", s, pkg.globals)    
}

#' @title Gets output directory
#' 
#' @description Returns the value of module variable for 
#' storing the current output directory. 
#' 
#' @export
getOutputDir<- function() {
  return(get("pkg.outputdir", pkg.globals))
}

#' @title getLogDir()
#' 
#' @description Returns the value for log directory
#' 
#' @export
getLogDir<- function() {
  return(paste0(getOutputDir(),"Log/"))
}

#' @title Create output directory
#' 
#' @description A simple function to make a directory to save the 
#' model's data.
#' 
#' @details Create the, if required, the directory to save the 
#' output data generate by the model. It is intended for internal 
#' use.
#' 
#' @export 
createOutputDir<- function() {
  lambda<- function(d) {
    if(!dir.exists(d)) {
      dir.create(d)
    }
  }
  
  ## -- Create required directories
  lambda(getOutputDir())
  lambda(getLogDir())
}



# Set the directory where repast model is installed
setBaseDir<- function(s) {
  assign("pkg.basedir", s, pkg.globals)  
}

# Gets the directory where repast model is installed
getBaseDir<- function() {
  return(get("pkg.basedir", pkg.globals))  
}

# Sets the directory where repast model is installed which normally
# is a subdirectory below installation base directory
setModelDir<- function(s) {
  assign("pkg.modeldir", s, pkg.globals)  
}

# Sets the directory where repast model is installed which normally
# is a subdirectory below installation base directory
getModelDir<- function() {
  return(get("pkg.modeldir", pkg.globals))  
}

# Sets the model's scenario directory
setScenarioDir<- function(s) {
  assign("pkg.scenariodir", s, pkg.globals)  
}

# Gets the model's scenario directory
getScenarioDir<- function() {
  return(get("pkg.scenariodir", pkg.globals))  
}

# Sets the model's lib directory
setModelLibDir<- function(s) {
  assign("pkg.modellibdir", s, pkg.globals)  
}

# Gets the model's lib directory
getModelLibDir<- function() {
  return(get("pkg.modellibdir", pkg.globals))  
}

# Traverse the lib dir to build up the classpath
repastlibs<- function() {
  libdir<- simphonylib()
  for(d in dir(libdir)) {
    # On bin dir we expect unpackaged class files
    bin<- paste0(libdir,paste0(d,"/bin"))
    lib<- paste0(libdir,paste0(d,"/lib"))
    # adding the bin dir to rjava classpath
    .jaddClassPath(bin)
    
    repastjars(paste0(libdir,d))
    repastjars(lib)
    repastjars(getModelLibDir())
  }
}

# Search for jar files inside lib dir and then add it to classpath
repastjars<- function(lib) {
  jars<- dir(lib,pattern="*.jar")
  for(j in jars) {
    jar<- paste0(paste0(lib,"/"),j)
    # adding jar file to classpath
    .jaddClassPath(jar)
  }
}


# ----- Exposed package API functions


#' @title jvm.set_parameters
#' 
#' @description Configures the jvm parameters
#' 
#' @details Set the underlying parameters for java virtual machine. The default 
#' values are "-server -Xms512m -Xmx1024m". These defaults can be changed 
#' to fit the model requirements.
#' 
#' @param s The paramter string to be passed to the underlying JVM 
#' 
#' @examples \dontrun{
#'    jvm.set_parameters("-server -Xms512m -Xmx2048m")}
#'    
#' @export 
jvm.set_parameters<- function(s) {
  assign("pkg.java.parameters", s, pkg.globals)  
}

#' @title jvm.get_parameters
#' 
#' @description Returns the current java virtual machine parameters
#' 
#' @return A string with JVM parameters. 
#' 
#' @export
jvm.get_parameters<- function() {
  return(get("pkg.java.parameters", pkg.globals))  
}

#' @title Init R/JVM environment
#' 
#' @description Initialize rJava and repast environment with classpath. This function
#' is called internally and it is not meant to be used directlly.
#' 
#' @details The default parameters can be changed as needed calling the 
#' primitive \code{\link{jvm.set_parameters}} befor instantiating the model 
#' engine.
#' 
#' @examples \dontrun{
#'      jvm.init()}
#' 
#' @references
#' [1] rJava: Low-Level R to Java Interface. Low-level interface to Java VM 
#' very much like .C/.Call and friends. Allows creation of objects, 
#' calling methods and accessing fields.
#' 
#' @import rJava
jvm.init<- function() {
  # The default parameters can be changed as needed
  .jinit(parameters= jvm.get_parameters() )
  .jaddClassPath(enginejar())
  .jaddClassPath(paste0(getModelDir(),"/bin"))
  # ----- Repast base libraries
  repastlibs()
}

#' @title jvm.setOut
#' 
#' @description Set the System.out filed to a file
#' 
#' @param f The output file name
#' 
#' @examples \dontrun{
#'    jvm.setOut("/tmp/SysteOut.log")}
#'    
#' @import rJava
#' 
#' @export 
jvm.setOut<- function(f) {
  ## -- Create the output dir if required
  createOutputDir()
  
  my.f<- paste0(getLogDir(),f)
  ## -- Java calls to redirect the output to a file
  obj.fos<- new(J("java.io.FileOutputStream"),my.f)
  obj.ps<- new(J("java.io.PrintStream"),obj.fos)
  .jcall("java/lang/System","V", "setOut",obj.ps)    
}

#' @title jvm.resetOut
#' 
#' @description Reset the System.out filed value to console output
#' 
#' @examples \dontrun{
#'    jvm.resetOut()}
#'    
#' @import rJava
#' 
#' @export 
jvm.resetOut<- function() {
  obj.os<- .jfield("java/io/FileDescriptor", name="out")
  obj.fos<- new(J("java.io.FileOutputStream"),obj.os)
  obj.ps<- new(J("java.io.PrintStream"),obj.fos)
  .jcall("java/lang/System","V", "setOut",obj.ps)    
}

# ----- Wrapper functions for Engine class method calls


#' @title Engine
#' 
#' @description Creates an instance of Engine
#' 
#' @details This function creates an instance of Repast model wrapper 
#' class. Before invoking the function Engine, make sure that 
#' environment was correctly initialized.
#' 
#' @return An onject instance of Engine class
#' 
#' @export
Engine<- function() {
  return(new(J(engineclazz())))
}

#' @title Engine.LoadModel
#' 
#' @description Loads the model's scenario files
#' 
#' @details This function loads the scenario of a Repast Model and 
#' initialize de model.
#' 
#' @param e An engine object instance
#' @param f The full path of scenario directory 
#' 
#' @export 
Engine.LoadModel<- function(e,f) {
  .jcall(e,"V", "LoadModel",f)  
}

#' @title Engine.SetAggregateDataSet
#' 
#' @description Sets the model's dataset
#' 
#' @details Configure a dataset with the desired output values 
#' to be "drained" by the function Engine.GetModelOutput. 
#' 
#' @param e An engine object instance
#' @param k The repast model's data set name
#' 
#' @examples \dontrun{
#'    d<- "C:/usr/models/your-model-directory"
#'    m<- Model(d)
#'    setAggregateDataSet(m,"dataset-name")}
#'    
#' @export
Engine.SetAggregateDataSet<- function(e,k) {
  .jcall(e,"V","ModelDataSet",k)
}

#' @title Engine.getParameterNames
#' 
#' @description Get the parameter names
#' 
#' @details Returns the names of all declared model's parameters in
#' the parameter.xml file in the scenario directory.
#' 
#' @param e An engine object instance
#' 
#' @return A collection of parameter names
#' 
#' @export
Engine.getParameterNames<- function(e) {
  names<- .jcall(e,"[S","getParameterNames")
  return(names)
} 

#' @title Engine.getParameter
#' 
#' @description The function gets the value of model 
#' parameter \code{k} as java.lang.Object
#' 
#' @param e An engine object instance
#' @param k The parameter name
#' 
#' @return The parameter value
#' 
#' @export
Engine.getParameter<- function(e,k) {
  v<- .jcall(e,"Ljava/lang/Object;","getParameter",k)
  return(v)
}

#' @title Engine.getParameterType
#' 
#' @description Returns the declared type of a Repast 
#' model parameter
#' 
#' @param e An engine object instance
#' @param k The parameter name
#' 
#' @return The parameter type string
#' 
#' @export 
Engine.getParameterType<- function(e,k) {
  v<- .jcall(e,"Ljava/lang/String;","getParameterType",k)
}

#' @title Engine.getParameterAsString 
#' 
#' @description Get the value of model parameter \code{k} 
#' as \code{java.lang.String}
#' 
#' @param e An engine object instance
#' @param k The parameter name
#' 
#' @return The parameter value as string
#' 
#' @export
Engine.getParameterAsString<- function(e,k) {
  v<- .jcall(e,"Ljava/lang/String;","getParameterAsString",k)
  return(v)
}

#' @title Engine.getParameterAsNumber
#' 
#' @description Get the value of model parameter 
#' \code{k} as \code{java.lang.Number}
#' 
#' @param e An engine object instance
#' @param k The parameter name
#' 
#' @return The parmeter value as number
#' 
#' @export
Engine.getParameterAsNumber<- function(e,k) {
  v<- .jcall(e,"Ljava/lang/Number;","getParameterAsNumber",k)
  return(v)
}

#' @title Engine.getParameterAsDouble
#' 
#' @description Get the value of model parameter \code{k} as \code{java.lang.Double}
#' 
#' @param e An engine object instance
#' @param k The parameter name
#' 
#' @return The parmeter value as double
#' 
#' @export
Engine.getParameterAsDouble<- function(e,k) {
  v<- .jcall(e,"D","getParameterAsDouble",k)
  return(v)
}

#' @title Engine.setParameter
#' 
#' @description Set the value of model parameter 
#' 
#' @param e An engine object instance
#' @param k The parameter name
#' @param v The parameter value
#' 
#' @export
Engine.setParameter<- function(e,k,v) {
  # Map the R type system to java object
  switch(Engine.getParameterType(e,k),
         java.lang.String = {
           value<- new(J("java.lang.String"), as.character(v))  
         },
         
         double = {
          value<- new(J("java.lang.Double"), as.double(v))
         },
         
         int = {
           value<- new(J("java.lang.Integer"), as.integer(v))
         },
         
         boolean = {
           value<- new(J("java.lang.Boolean"), as.logical(v))
         })
  # Invoke the setParamter method
  .jcall(e,"V","setParameter",k,value)
}

#' @title Engine.endAt
#' 
#' @description Configure the maximun simulated time for 
#' the current model run
#' 
#' @param e An engine object instance
#' @param v The number of Repast time ticks 
#' 
#' @export
Engine.endAt<- function(e,v) {
  .jcall(e,"V","endAt",v)
}

#' @title Returns the model id
#' 
#' @description This function provides a wrapper to the method getId() 
#' from repast context. The id is basically a String with the currently 
#' instantiated model name.
#' 
#' @param e An engine object instance
#' 
#' @export
Engine.getId<- function(e) {
  id<- .jcall(e,"S","getId")
  if(nchar(id) == 0) stop("Model not initilized.")
  return(id)
}

#' @title Engine.RunModel
#' 
#' @description Performs the execution of Repast model
#' 
#' @param e An engine object instance
#' 
#' @export
Engine.RunModel<- function(e) {
  .jcall(e,"V","RunModel")
}

#' @title Engine.GetModelOutput
#' 
#' @description Gets the model output data as a CSV String array. 
#' Calls the engine method GetModelOutput to drain model output
#' data.
#' 
#' @param e An engine object instance
#' 
#' @return An array of strings containing the model's output
#' 
#' @examples \dontrun{
#'    d<- "c:/usr/models/your-model-directory"
#'    m<- Model(d)
#'    csv<- Engine.GetModelOutput(m)}
#' @importFrom utils read.csv
#' @export
Engine.GetModelOutput<- function(e) {
  .jcall(e,"[S","GetModelOutput")
}

#' @title Engine.Finish
#' 
#' @description Performs a cleanup on a engine instance.Finalize 
#' and destroy repast controller data.
#' 
#' @param e An engine object instance
#' 
#' @export
Engine.Finish<- function(e) {
  .jcall(e,"V","cleanUpBatch")
}

#' @title Set the log level to INFO
#' 
#' @description Configures the underlying logging system
#' 
#' @export
Logger.setLevelInfo<- function() {
  logger<- J("org.haldane.rrepast.RepastEngineLogger")
  .jcall(logger,"V","setLevelInfo")
}

#' @title Set the log level to WARNING
#' 
#' @description Configures the underlying logging system
#' 
#' @export
Logger.setLevelWarning<- function() {
  logger<- J("org.haldane.rrepast.RepastEngineLogger")
  .jcall(logger,"V","setLevelWarning")
}

#' @title ShowModelPaths
#' 
#' @description Prints the paths. Shows the directories 
#' currently used to load model scenario and lib. The output of 
#' this function is informational only and can be used to check 
#' whether model data is being loaded properly from 
#' correct locations.
#' 
#' @examples \dontrun{
#'    ShowModelPaths()}
#'    
#' @export
ShowModelPaths<- function() {
  print(paste("Install dir.= ",getBaseDir()))
  print(paste("Model dir...= ",getModelDir()))
  print(paste("Scenario dir= ",getScenarioDir()))
  print(paste("Model lib...= ",getModelLibDir()))
}

#' @title ShowClassPath
#' 
#' @description Shows the current classpath
#' 
#' @return the current setting of JVM classpath
#' 
#' @examples \dontrun{
#'    ShowClassPath()}
#'    
#' @export
ShowClassPath<- function() {
  .jclassPath()
}

### --- Progress Bar functions

#' @title PB.set
#' 
#' @description  Ses the progress bar descriptor
#' 
#' @param obj -- The progress bar descriptor
#' 
#' @export
PB.set<- function(obj) {
  assign("pkg.progressbar", obj, pkg.globals)      
}

#' @title PB.get
#' 
#' @description  Gets the the progress bar descriptor
#' 
#' @export
PB.get<- function() {
  return(get("pkg.progressbar", pkg.globals))
}

#' @title PB.enable
#' 
#' @description  Enables the progress bar visualization
#' 
#' @export
PB.enable<- function() {
  assign("pkg.progressbar.enabled", TRUE, pkg.globals)      
}

#' @title PB.disable
#' 
#' @description  Disable the progress bar visualization
#' 
#' @export
PB.disable<- function() {
  assign("pkg.progressbar.enabled", FALSE, pkg.globals)      
}

#' @title PB.isEnabled
#' 
#' @description Returns the global value indicating if progress bar 
#' is enabled.
#' 
#' @return Boolean TRUE if progress bar must be shown
#' 
#' @export
PB.isEnabled<- function() {
  return(get("pkg.progressbar.enabled", pkg.globals))
}

#' @title PB.init
#' 
#' @description Initialize progress bar for model
#' execution.
#' 
#' @param psets -- The total number of paramter sets being simulated
#' @param replications -- The number of replications per simulation round
#' 
#' @importFrom utils setTxtProgressBar txtProgressBar
#' 
#' @export
PB.init<- function(psets, replications) {
  
  ## -- Check if init function has already been called from RunExperiment
  if(length(grep("(RunExperiment\\s*\\()|(Run\\s*\\()",sys.calls())) == 2) {
    ##print(grep("(RunExperiment\\s*\\()|(Run\\s*\\()",sys.calls(),value=TRUE))  
    return()
  }
  
  ## -- print(grep("(RunExperiment\\s*\\()|(Run\\s*\\()",sys.calls(),value=TRUE))  
  ## -- print("PB.init")
  
  if(PB.isEnabled()) {
    total<- psets * replications
    pbar<- txtProgressBar(min = 0, max = total, style = 3)
    pbar$pset<- 1
    pbar$replications<- replications
    PB.set(pbar)  
  }
  return(sys.calls())
}

#' @title PB.close
#' 
#' @description Close the progress bar descriptor
#' 
#' @export
PB.close<- function() {
  ## -- Check if init function has already been called from RunExperiment
  if(length(grep("(RunExperiment\\s*\\()|(Run\\s*\\()",sys.calls())) == 2) {
    return()
  }
  
  if(PB.isEnabled()) {
    pbar<- PB.get()
    if(!is.null(pbar)) {
      close(pbar)
      PB.set(NULL)  
    } else {
      stop("Progress bar has not been initialized!")
    }
  }
}

#' @title PB.pset
#' 
#' @description Update pset value
#' 
#' @param v The current parameter set being simulated
#' 
#' @export
PB.pset<- function(v) {
  if(PB.isEnabled()) {
    pbar<- PB.get()
    if(!is.null(pbar)) {
      pbar$pset<- v
      PB.set(pbar)  
    }  else {
      stop("Progress bar has not been initialized!")
    }
  }
}

#' @title PB.update
#' 
#' @description Update progress bar
#' 
#' @param r The current replication number
#' 
#' @importFrom utils setTxtProgressBar txtProgressBar
#' 
#' @export
PB.update<- function(r) {
  if(PB.isEnabled()) {
    pbar<- PB.get()
    if(!is.null(pbar)) {
      setTxtProgressBar(pbar, (pbar$pset-1)*pbar$replications + r)

    }  else {
      stop("Progress bar has not been initialized!")
    }
  }
}

#' @title The easy API for model initilization
#' 
#' @description Instantiate a repast model from the model dir without
#' loading the scenario file.
#' 
#' @details This is the entry point for model execution. Typically 
#' any model execution will start with this function which encapsulates 
#' all low level calls for model initialization. In order to perform 
#' simulations with repast from R code only \code{Model} and a 
#' few more function calls are required: \code{\link{Load}}, 
#' \code{\link{Run}}. Finally the output of model is managed with 
#' functions \code{\link{GetResults}} and \code{\link{SaveSimulationData}}.
#' 
#' @param modeldir The installation directory of some repast model
#' @param maxtime The total simulated time
#' @param dataset The name of any model aggregate dataset
#' @param load If true instantiate model and load scenario
#' 
#' @return Returns the instance of repast model
#' @examples \dontrun{
#'    d<- "C:/usr/models/your-model-directory"
#'    m<- Model(d)}
#'    
#' @references
#' [1] North, M.J., N.T. Collier, and J.R. Vos, "Experiences Creating Three Implementations of the Repast Agent Modeling Toolkit," ACM Transactions
#' on Modeling and Computer Simulation, Vol. 16, Issue 1, pp. 1-25, ACM,
#' New York, New York, USA (January 2006).
#' @export
Model<- function(modeldir="",maxtime=300,dataset="none", load=FALSE) {
  if(dir.exists(modeldir)) {
    # Configure all required directory based on their default locations
    configModelDirs(modeldir)
    
    # Initilialized JVM
    # Setup classpath inferring default values from modeldir
    jvm.init()
    
    # Creating an engine instance
    e<- Engine()
    
    # Configure the total amount of time to be simulated
    Engine.endAt(e,maxtime)
    
    # Configure the dataset
    Engine.SetAggregateDataSet(e,dataset)
    
    if(load == TRUE) {
      Load(e)
    }
    
    return(e)
  } else {
    stop(paste0("The model directory does not exist: ", modeldir))
  }
}

#' @title The Scenario loader
#' 
#' @description Loads the model's scenario. This function must be 
#' called before running the model.
#' 
#' @examples \dontrun{
#'    d<- "C:/usr/models/your-model-directory"
#'    m<- Model(d)
#'    Load(m)}
#' 
#' @param e An engine object instance
#' 
#' @export
Load<- function(e) {
  Engine.LoadModel(e,getScenarioDir())
  setId(Engine.getId(e))
}

#' @title Run simulations
#' 
#' @description This function executes the time steps of an 
#' instantiated model. The number of replications of model 
#' runs can be specified by the function parameter. The seed 
#' parameter may be omitted and will be generated internally. 
#' If provided, the seed collection, must contain the same 
#' number of \code{r} parameter. 
#'
#' @param e An engine object instance
#' @param r The number of experiment replications
#' @param seed The random seed collection
#' 
#' @return The model output dataset
#' 
#' @examples \dontrun{
#'    d<- "C:/usr/models/your-model-directory"
#'    m<- Model(d)
#'    Load(m)
#'    Run(m,r=2) # or Run(m,r=2,seed=c(1,2))}
#'    
#' @importFrom stats runif
#' @export
Run<- function(e,r=1,seed=c()) {
  # The default behaviour is if seed set was
  # not provided generate a suitable set of 
  # random seeds for the number of replications.
  if(length(seed) == 0) {
    seed= runif(r,-10^8,10^8)  
  } else if(length(seed) != r) {
    stop("The provided set of random numbers doesn't match replications!")
  }
  
  # Gets the current set of parameters
  p<- GetSimulationParameters(e)
  
  # Clear result repository
  ClearResults()
  
  SetResultsParameters(p)
  
  ## --- Progress Bar initialization
  PB.init(1, r)
  
  results<- c()
  
  for(i in 1:r) {
    ## --- Setting the random seed for experiment replication
    Engine.setParameter(e,getKeyRandom(),as.integer(seed[i]))

    ## --- Pass the control to Repast to run simulation      
    Engine.RunModel(e)

    data<- GetOutput(e)
    
    ## --- Just for debug: print(paste0("run= ", i, " / rows=", nrow(data)))
    
    # Sets the current run number
    data$run<- i
    
    # Add replication data 
    AddResults(data)
    
    results<- rbind(results,data)
    
    # Update progress bar
    PB.update(i)
  }
  ## --- Progress Bar clean up
  PB.close()
  
  return(results)
}

#' @title Run an experimental setup
#' 
#' @description Run the model multiple times for different parameters
#' given by design matrix function parameter.
#' 
#' @details The FUN function must return zero for perfect fit and values 
#' greater than zero otherwise.
#'
#' @param e An engine object instance
#' @param r The number of experiment replications
#' @param design The desing matrix holding parameter sampling
#' @param FUN THe calibration function.
#' 
#' @examples \dontrun{
#'    my.cost<- function(params, results) { # your best fit calculation, being 0 the best metric.  }
#'    d<- "c:/usr/models/your-model-directory"
#'    m<- Model(d,dataset="ds::Output")
#'    Load(m)
#'    f<- AddFactor(name="cyclePoint",min=40,max=90)
#'    f<- AddFactor(factors=f, name="conjugationCost",min=1,max=80)
#'    d<- LatinHypercube(factors=f)
#'    p<- GetSimulationParameters(e)
#'    exp.design<- BuildParameterSet(d,p)
#'    v<- RunExperiment(e,r=1,exp.design,my.cost) }
#'    
#' @return A list with output and dataset
#' 
#' @export
RunExperiment<- function(e, r=1, design, FUN) {
  paramset<- c()
  output<- c()
  dataset<- c()
  
  psets<- nrow(design)
  
  ## --- Progress Bar initialization
  PB.init(psets, r)
  
  for(i in 1:psets) {
    d<- design[i,]
    
    # -- Set parameters for next model 'Run'
    SetSimulationParameters(e, d)
    
    # -- Update progress bar pset value
    PB.pset(i)
    
    # -- Run model with current parameter set
    Run(e,r)
    
    results<- GetResults()
    
    # -- The user provided calibration function.
    # -- Calibration function must return 0 for perfect fit between 
    # -- observed and experimental data.
    calibration<- FUN(d,results)
    
    if(is.null(calibration)) {
      stop("Invalid user provided calibration function!")  
    }
    
    # -- The 'pset' is the paramter set id
    pset<- i
    
    paramset<- rbind(paramset,cbind(pset,d))
    output<- rbind(output,cbind(pset,calibration))
    dataset<- rbind(dataset,cbind(pset,results))
  }
  
  ## --- Progress Bar clean up
  PB.close()
  
  return(list(paramset=paramset, output=output, dataset=dataset))
}

#' @title Helper function to get experiment \code{paramset}
#' 
#' @description The RunExperiment function returns a list holding
#' the \code{paramset}, \code{output} and \code{dataset} collection.
#' The \code{paramset} collection contains the parameters used for 
#' running the experimental setup. The \code{output} has the results 
#' from user provided calibration function. The \code{dataset} 
#' collection has the raw output of 'Repast' aggregated dataset.
#' 
#' @examples \dontrun{
#'    d<- "C:/usr/models/your-model-directory"
#'    m<- Model(d)
#'    ...
#'    e<- RunExperiment(e,r=1,exp.design,my.cost)
#'    p<- getExperimentParamSet(e)}
#' 
#' @param e The experiement object returned by \code{\link{RunExperiment}} 
#' 
#' @return The reference to \code{output} container.
#' @export
getExperimentParamSet<- function(e) {
  v<- e$paramset
  return(v)
}

#' @title Helper function to get experiment \code{output}
#' 
#' @description The RunExperiment function returns a list holding
#' the \code{paramset}, \code{output} and \code{dataset} collection.
#' The \code{paramset} collection contains the parameters used for 
#' running the experimental setup. The \code{output} has the results 
#' from user provided calibration function. The \code{dataset} 
#' collection has the raw output of 'Repast' aggregated dataset.
#' 
#' @param e The experiement object returned by \code{\link{RunExperiment}} 
#' 
#' @return The reference to \code{output} container.
#' @export
getExperimentOutput<- function(e) {
  v<- e$output
  return(v)
}

#' @title Helper function to get experiment \code{dataset}
#' 
#' @description The RunExperiment function returns a list holding
#' the \code{paramset}, \code{output} and \code{dataset} collection.
#' The \code{paramset} collection contains the parameters used for 
#' running the experimental setup. The \code{output} has the results 
#' from user provided calibration function. The \code{dataset} 
#' collection has the raw output of 'Repast' aggregated dataset.
#' 
#' @param e The experiement object returned by \code{\link{RunExperiment}} 
#' 
#' @return The reference to \code{dataset} container.
#' @export
getExperimentDataset<- function(e) {
  v<- e$dataset
  return(v)
}

#' @title Gets the output
#' 
#' @description  Returns the results of a model a data.frame from the last
#' RUN. Should be used only if model replication is equal to 1,
#' otherwise GetResults must be used.
#' 
#' @param e An engine object instance
#' 
#' @return Returns a data.frame with output data
#' 
#' @examples \dontrun{
#'    d<- "C:/usr/models/your-model-directory"
#'    m<- Model(d)
#'    ...
#'    data<- GetOutput(m)}
#'    
#' @importFrom utils read.csv
#' @export
GetOutput<- function(e) {
  c<- textConnection(Engine.GetModelOutput(e))
  read.csv(c)
}

#' @title Set parameters for running model
#' 
#' @description Modify the repast model parameters with 
#' values provided in parameter 'p' which is a data frame
#' with just one row.
#' 
#' @param e An engine object instance
#' @param p A data frame with simulation parameters
#' 
#' @export
SetSimulationParameters<- function(e, p) {
  if(is.null(e)) {
    stop("Engine object is null!")  
  }

  for(key in names(p)) {
    value<- p[1,key]
    if(is.factor(value)) {
      value<- levels(value)
    } 
    
    ## Modify the default set of parameters
    SetSimulationParameter(e, key, value)
  }
}

#' @title SetSimulationParameter 
#' 
#' @description Modify model's default parameter collection
#' 
#' @param e An engine object instance
#' @param key The paramter name
#' @param value The parameter value
#' 
#' @export
SetSimulationParameter<- function(e, key, value) {
  if(is.null(e)) {
    stop("Engine object is null!")  
  }
  
  keys<- names(GetSimulationParameters(e))
  
  # Verify that "key" is a valid model parameter
  if(key %in% keys) {
  # Try to coerce the value to a type for safety
    switch(typeof(value),
      double = {
        #print(paste0("double", key,"<- ",value))
        value<- as.double(value)
      },
      
      integer = { 
        #print(paste0("integer", key,"<- ",value))
        value<- as.integer(value)
      },
      
      character = { 
         #print(paste0("character", key,"<- ",value))
         value<- as.character(value)
       })
    Engine.setParameter(e,key,value)          
  }
}

#' @title Gets the simulation parameters
#' 
#' @description Returns a dataframe with the current set of input 
#' parameters for the last model run.
#' 
#' @param e An engine object instance
#' 
#' @return A data frame with simulation parameters
#'
#' @export
GetSimulationParameters<- function(e) {
  keys<- ""
  values<- ""
  names<- Engine.getParameterNames(e)
  for(n in names) {
    v<- Engine.getParameterAsString(e,n)
    if(nchar(keys) == 0){
      keys<- n
      values<- v
    } else {
      keys<- paste0(keys,",",n)
      values<- paste0(values,",",v)
    }
  }
  b<- rbind(keys,values)
  c<- textConnection(b)
  read.csv(c)
} 

#' @title Clear the results data.frame
#' 
#' @description This function is called automatically every
#' time Run method is called.
#'
#' @export
ClearResults<- function() {
  assign("pkg.results", data.frame(), pkg.globals)   
  assign("pkg.parameters", data.frame(), pkg.globals)   
}

#' Returns the model results
#'
#' @export
GetResults<- function() {
  return(get("pkg.results", pkg.globals))
}

#' Stores a data.frame 
#' 
#' @param d A data frame containing one replication data
#'
#' @export
SetResults<- function(d) {
  assign("pkg.results", d, pkg.globals)     
}

#' @title Concatenate results of multiple runs
#' 
#' @description This function stores the output 
#' of the last model execution and it is intended 
#' to be used internally.
#' 
#' @param d A data frame containing one replication data
#'
#' @export
AddResults<- function(d) {
  r<- GetResults()
  SetResults(rbind(r,d))
}

#' @title Gets the parameters
#' 
#' @description Returns the current set of paramters used 
#' for the last model run.
#'
#' @return A data.frame with parameters of the model.
#'
#' @export
GetResultsParameters<- function() {
  return(get("pkg.parameters", pkg.globals))
}

#' @title Sets the parameters
#' 
#' @description Save the current set of paramters used 
#' for the last model run.
#' 
#' @param d A data.frame with parameter values
#'
#' @export
SetResultsParameters<- function(d) {
  assign("pkg.parameters", d, pkg.globals)     
}

#' @title Saving simulation output
#' 
#' @description Saves the simulation results of last call to Run(e)
#' function.
#' 
#' @details The model must have been initialized or user must call 
#' \code{setId} explicitelly.
#' 
#' @param as The desired output type, must be csv or xls
#' @param experiment The experiment output
#' 
#' @return The id of saved data
#' 
#' @importFrom xlsx write.xlsx
#' @importFrom digest digest
#' @importFrom utils write.csv
#' @export
SaveSimulationData<- function(as="csv", experiment=NULL) {
  # Creating output dir if needed
  createOutputDir()
  filename<- getId()
  if(is.na(filename)) {
    stop("Model was not initialized correctly!")
  }
  
  paramset<- NULL
  output<- NULL
  dataset<- NULL
  
  if(!is.null(experiment)) {
    paramset<- getExperimentParamSet(experiment)
    output<- getExperimentOutput(experiment)
    dataset<- getExperimentDataset(experiment)
    
  } else {
    # The parameters of current simultation output
    paramset<- GetResultsParameters()
    
    # The results of simulation run
    dataset<- GetResults()  
  }
  
  
  
  hash <- digest(Sys.time(), algo="crc32")
  f0<- paste0(getOutputDir(),tolower(filename),"-paramset-",hash)
  f1<- paste0(getOutputDir(),tolower(filename),"-output-",hash)
  f2<- paste0(getOutputDir(),tolower(filename),"-dataset-",hash)
  
  switch(as,
         csv = {
           f0<- paste0(f0,".csv")
           f1<- paste0(f1,".csv")
           f2<- paste0(f2,".csv")
           write.csv(paramset, f0, row.names=FALSE)
           if(!is.null(output)) {
             write.csv(output, f1, row.names=FALSE)  
           }
           write.csv(dataset, f2, row.names=FALSE)
         },
         xls = { 
           f0<- paste0(f0,".xlsx")
           f1<- paste0(f1,".xlsx")
           f2<- paste0(f2,".xlsx")
           
           write.xlsx(paramset, f0)
           if(!is.null(output)) {
             write.xlsx(output, f1)
           }
           write.xlsx(dataset, f0)
         })
  return(hash)
}

##
## ----- Below sensitivity analysis functions
##

#' @title Adds a paramter to factor collection
#' 
#' @description Builds up the factor collection.
#' 
#' @param factors The current factor collection
#' @param lambda The function to apply FUN(p,min,max)
#' @param name The name of factor
#' @param min The minimun of parameter p
#' @param max The maximun of parameter p
#' 
#' @examples \dontrun{
#'    f<- AddFactor(name="Age",min=20,max=60)
#'    f<- AddFactor(factors=f, name="Weight",min=50,max=120)}
#' 
#' @return The collection of created factors
#'
#' @export
AddFactor<- function(factors=c(), lambda="qunif",name, min, max) {
  if(max < min) {
    stop("Invalid factor range!")
  }
  
  # if parameter already existe replace the current value
  rrow<- c(lambda=lambda,name=name,min=min,max=max)
  rownames(rrow)<- NULL
  if(length(factors) > 0 && factors[,"name"] == name) {
    i<- which(factors[,"name"] == name)  
    factors[i,]<- c(rrow)
  } else {
    factors<- rbind(factors,c(rrow))
  }
  return(factors)
}

#' @title Get the number of factors
#' 
#' @description Returns the total number of factors
#' 
#' @param factors A collection of factors created with AddFactor
#' 
#' @return The number of parameters in factors collection
#'
#' @export
GetFactorsSize<- function(factors) {
  n<- nrow(factors)
  if(is.null(n)) n<- 0
  return(n)
}

#' @title Corrects the LHS design matrix
#' 
#' @description Correct the LHS sampling matrix for a 
#' specific range applying the lambda function. The default
#' value of 'lambda' is 'qunif'.
#' 
#' @param design The LHS design matrix
#' @param factors THe collection of factors
#' 
#' @return The corrected design matrix
#'
#' @export
ApplyFactorRange<- function(design, factors) {
  k<- GetFactorsSize(factors)
  d<- sapply(1:k, function(p) {match.fun(factors[p,"lambda"])(design[,p],as.numeric(factors[p,"min"]),as.numeric(factors[p,"max"]))})
  
  if(is.null(nrow(d))) {
    ## --- Handle the case where sample size is 1
    d<- as.data.frame(t(d))  
  } else {
    ## --- Handle the case where sample size > 1
    d<- as.data.frame(d)  
  }

  names(d)<- factors[,"name"]
  return(d)
}

#' @title Builds the simulation parameter set
#' 
#' @description Merges the design matrix with parameters which 
#' will be keep fixed along simulation runs.
#' 
#' @param design The experimental desing matrix for at least one factor
#' @param parameters All parameters of the repast model.
#' 
#' @return A data frame holding all parameters required for running the model
#' 
#' @examples \dontrun{
#'    modeldir<- "c:/usr/models/BactoSim(HaldaneEngine-1.0)"
#'    e<- Model(modeldir=modeldir,dataset="ds::Output")
#'    Load(e)
#'    
#'    f<- AddFactor(name="cyclePoint",min=40,max=90)
#     f<- AddFactor(factors=f, name="conjugationCost",min=1,max=80)
#'    
#'    p<- GetSimulationParameters(e)
#'    
#'    d<- LatinHypercube(factors=f)
#'    
#'    p1<- BuildParameterSet(d,p)}
#' 
#' @export
BuildParameterSet<- function(design, parameters) {
  v<- as.data.frame(design)
  tmp.p<- parameters
  for(n in names(v)) {
    # Drop parameters columns which are in design matrix
    tmp.p[n]<- NULL
  }

  # Now join two data frames
  for(i in 1:length(names(tmp.p))) {
    v<- cbind(tmp.p[i],v)      
  }
  return(v)
}

#' @title SequenceItem
#' 
#' @description Generate a sequence from min to max using an increment
#' based on the number of of elements in v
#' 
#' @param v A column of n x k design matrix
#' @param min The lower boundary of range
#' @param max The uper boundary of range
#' 
#' @return A sequence between min and max value
#' 
#' @export
SequenceItem<- function(v,min,max) {
  n<- length(v)
  delta<- (max-min)/(n-1)
  return(seq(min,max,delta))
}

#' @title df2matrix
#' 
#' @description This function converts data frames to matrix data type.
#' 
#' @param d The data frame
#' @param n The column names to be converted. Null for all data frame columns
#' 
#' @return The data frame converted to a matrix
#' 
#' @export
df2matrix<- function(d,n=c()) {
  if(length(n) == 0) {
    n<- names(d)
  }
  m<- c()
  for(k in n) {
    m<- cbind(m,as.matrix(d[,k]))
  }
  colnames(m)<- n
  return(m)
}

#' @title dfsumcol
#' 
#' @description Sum data frame columns but tho
#' 
#' @param d The data frame 
#' @param lst Skip columns included. Sum columns NOT included 
#' @param invert Sum only the columns included in \code{lst}
#' 
#' @return The original data frame with a new column (sum) holding the sum 
#' 
#' @export
dfsumcol<- function(d,lst=c(),invert=FALSE) {
  v<- as.data.frame(d)
  s<- NULL
  
  op<- "!"
  if(invert) {
    identity<- function(x) {x}
    op<- "identity"
  }

  for(key in colnames(v)) {
    if(match.fun(FUN=op)(toupper(key) %in% toupper(lst))) {
      if(is.null(s)){
        s<- v[, key]
      } else {
        s<- s + v[, key]
      }
    }
  }  
  v$total<- s
  
  ## Return the same type of original data
  if(is.matrix(d)) {
    v<- as.matrix(v)  
  }
  
  return(v)
}

#' @title dffilterby
#' 
#' @description Selects a subset of a data frame, filtering by 
#' column values.
#' 
#' @param d The data frame holding data to be filtered
#' @param key The column name for selection valuas
#' @param values The collection of values used to filter the data set
#' 
#' @return The filtered data set
#' 
#' @export
dffilterby<- function(d, key, values=c()) {
  d<- as.data.frame(d)
  
  o<- c()
  for(v in values) {
    o<- rbind(o,d[d[,colnames(d) == key] %in% v,])    
  }
  return(o)
}

#' @title dfround
#' 
#' @description Round all numeric columns of a data frame
#' 
#' @param d The data frame 
#' @param p The number of decimal digits to be keept
#' 
#' @return A data frame with rounded columns
#' 
#' @export
dfround<- function(d, p) {
  return(sapply(d[,sapply(d,is.numeric)],round,digits=p))
}

# ----- DoE - Design of Experiments  
# ----- AoE - Analysis of Experimental Data

#' @title AoE.RMSD
#' 
#' @description  A simple Root-Mean-Square Deviation 
#' calculation.
#' 
#' @param xs The simulated data set
#' @param xe The experimental data set
#' 
#' @return The RMSD value for provided datasets
#' @export
AoE.RMSD<- function(xs, xe) {
  return(sqrt(mean((xs - xe)^2, na.rm = TRUE)))
}

#' @title AoE.MAE
#' 
#' @description Calculates the average-error
#' magnitude (MAE)
#' 
#' @param xs The simulated data set
#' @param xe The experimental data set
#' 
#' @return The MAE value for provided datasets
#' @export
AoE.MAE<- function(xs, xe) {
  return(mean((abs(xs - xe)), na.rm = TRUE))
}

#' @title AoE.NRMSD
#' 
#' @description  A simple Normalized Root-Mean-Square 
#' Deviation calculation using max and min values.
#' NRMSD = RMSD(x) / (max(x) - min(x))
#' 
#' @param xs The simulated data set
#' @param xe The experimental data set
#' 
#' @return The NRRMSD value for provided datasets
#' 
#' @export
AoE.NRMSD<- function(xs, xe) {
  ##divisor<- abs(max(xe,na.rm = TRUE) - min(xe,na.rm = TRUE))
  return(sqrt(mean((xs - xe)^2, na.rm = TRUE))/mean(xe, na.rm = TRUE))
}

#' @title AoE.CoV
#' 
#' @description A simple funcion for calculate the 
#' Coefficient of Variation
#' 
#' @param d The data collection
#' @return The coefficient of variation for data
#' 
#' @importFrom stats sd
#' 
#' @export
AoE.CoV<- function(d) {
  return((sd(d,na.rm = TRUE)/mean(d,na.rm = TRUE)) * 100)
}

#' @title AoE.ColumnCoV
#' 
#' @description This function Calculates the relative squared 
#' deviation (RSD or CoV) for an used provided column name \code{key}
#' in the parameter \code{dataset}. 
#' 
#' @param dataset A model output dataset
#' @param key Column name from output dataset
#' 
#' @return A data frame with Coefficient of variations
#' 
#' @export
AoE.ColumnCoV<- function(dataset, key) {
  m.run<- dataset$run
  if(is.null(m.run)) {
    stop("The dataset is not an instance of model output!")
  }
  
  result<- c()
  m.max<- max(m.run)
  
  for(i in 1:m.max) {
    m.data<- with(dataset,dataset[run %in% seq(1,i), key])
    result<- rbind(result,cbind(i, AoE.CoV(m.data)))
  }
  result<- as.data.frame(result)
  names(result)<- c("sample","RSD")
  return(result)
}

#' @title AoE.Stability
#' 
#' @description This function verifies the stability 
#' of CoV for all columns given by parameter \code{keys}
#' or all dataset columns if keys is empty.
#' 
#' @param dataset A model output dataset
#' @param keys A list of column names
#' 
#' @return A data frame with Coefficient of variations
#' 
#' @export
AoE.Stability<- function(dataset, keys=c()) {
  if(length(keys) == 0) {
    keys<- setdiff(names(dataset), c("pset","random_seed","run","Time"))   
  }
  
  results<- c()
  for(k in keys) {
    v<- AoE.ColumnCoV(dataset,k)
    v$group<- k
    results<- rbind(results,v)
  }
  return(results)
}
  
#' @title AoE.Base
#' 
#' @description The Design Of Experiments Base function
#' 
#' @param m The base design matrix
#' @param factors A subset of model parameters
#' @param fun The function which will be applied to m
#'
#' @return The design matrix
#'
#' @export
AoE.Base<- function(m, factors=c(), fun=NULL) {
  k<- GetFactorsSize(factors)
  if(k == 0) {
    stop("Empty factor collection!")
  }
  
  tmp.factors<- factors
  if(!is.null(fun)) {
    tmp.factors[,"lambda"]<- fun  
  }
  
  # --- Apply the desired range
  design<- ApplyFactorRange(m, tmp.factors)
  return(design)
}

#' @title AoE.LatinHypercube 
#' 
#' @description Generate a LHS sample for model parameters
#' 
#' @details Generate the LHS sampling for evaluating 
#' the parameters of a model.
#' 
#' @param n The number of samples
#' @param factors The model's parameters which will be evaluated
#' 
#' @return The LHS design matrix for provided parameters
#' 
#' @examples \dontrun{
#'  f<- AddFactor(name="cyclePoint",min=40,max=90)
#'  f<- AddFactor(factors=f, name="conjugationCost",min=1,max=80)
#'  d<- DoE.LatinHypercube(2,f)}
#' 
#' @importFrom lhs randomLHS
#' @export
AoE.LatinHypercube<- function(n=10, factors=c()) {
  k<- GetFactorsSize(factors)
  
  # --- Generate design matrix
  design<- AoE.Base(randomLHS(n, k), factors)
  return(design)
}

#' @title AoE.FullFactorial design generator
#' 
#' @description Generate a Full Factorial sampling for evaluating 
#' the parameters of a model.
#' 
#' @param n The number of samples
#' @param factors The model's parameters which will be evaluated
#' 
#' @return The Full Factorial design matrix for provided parameters
#' 
#' @examples \dontrun{
#'  f<- AddFactor(name="cyclePoint",min=40,max=90)
#'  f<- AddFactor(factors=f, name="conjugationCost",min=1,max=80)
#'  d<- AoE.FullFactorial(2,f)}
#' 
#' @export
AoE.FullFactorial<- function(n=10, factors=c()) {
  k<- GetFactorsSize(factors)
  
  # --- calculate n for the aproximate number of samples
  n<- round(n^(1/k))
  
  # --- Generate design matrix
  design<- AoE.Base(matrix(nrow = n, ncol = k, seq(1,n)), factors, "SequenceItem")
  design<-  expand.grid(design)
  return(design)
}

#' @title AoE.RandomSampling experiment desing generator
#' 
#' @description Generate a Simple Random Sampling experiment design
#' matrix.
#' 
#' @param n The number of samples
#' @param factors The model's parameters which will be evaluated
#' 
#' @return The random sampling design matrix 
#' 
#' @examples \dontrun{
#'  f<- AddFactor(name="cyclePoint",min=40,max=90)
#'  f<- AddFactor(factors=f, name="conjugationCost",min=1,max=80)
#'  d<- AoE.RandomSampling(2,f)}
#' 
#' @export
AoE.RandomSampling<- function(n=10, factors=c()) {
  k<- GetFactorsSize(factors)
  m<- c()
  for(i in 1:k) {
    m<- cbind(m,runif(n))  
  }
  design<- AoE.Base(m, factors)
  return(design)
}

#' @title AoE.Morris 
#' 
#' @description This is a wrapper for performing Morris's  screening
#' method on repast models. We rely on morris method from sensitivity 
#' package.
#' 
#' @param k The factors for morris screening.
#' @param p The number of levels for the model's factors.
#' @param r Repetitions. The number of random sampling points of Morris Method.
#' 
#' @references Gilles Pujol, Bertrand Iooss, Alexandre Janon with contributions from Sebastien Da Veiga, Jana Fruth,
#' Laurent Gilquin, Joseph Guillaume, Loic Le Gratiet, Paul Lemaitre, Bernardo Ramos and Taieb Touati (2015).
#' sensitivity: Sensitivity Analysis. R package version 1.11.1.
#' https://CRAN.R-project.org/package=sensitivity
#' 
#' @importFrom sensitivity morris
#' @export
AoE.Morris<- function(k=c(),p=5,r=4) {
  
  k.v<- GetFactorsSize(k)
  if(k.v == 0) {
    stop("Empty factor collection!")
  }
  
  p.min<- as.numeric(k[,"min"])
  p.max<- as.numeric(k[,"max"])
  p.design<- list(type = "oat", levels = p, grid.jump = ceiling(p/2))  
  v<- morris(NULL, k[,"name"], r, p.design, p.min, p.max, scale=TRUE)
  return(v)
}

#' @title AoE.GetMorrisOutput 
#' 
#' @description  Returns a dataframe holding the Morris 
#' result set
#' 
#' @param obj A reference to a morris object instance
#' 
#' @return The results of Morris method
#' 
#' @importFrom stats sd
#' @export
AoE.GetMorrisOutput<- function(obj) {
  mu <- apply(obj$ee, 2, mean)
  mu.star <- apply(obj$ee, 2, function(x) mean(abs(x)))
  sigma <- apply(obj$ee, 2, sd) 
  m<- t(rbind(mu,mu.star,sigma))
  tmp<- as.data.frame(m,row.names=seq(1,nrow(m)))
  tmp$group<- rownames(m)
  return(tmp)
}

#' @title AoE.Sobol
#' 
#' @description This is a wrapper for performing Global Sensitivity
#' Analysis using the Sobol Method provided by sensitivity 
#' package.
#' 
#' @details This function is not intended to be used directly from 
#' user programs.
#' 
#' @references Gilles Pujol, Bertrand Iooss, Alexandre Janon with contributions from Sebastien Da Veiga, Jana Fruth,
#' Laurent Gilquin, Joseph Guillaume, Loic Le Gratiet, Paul Lemaitre, Bernardo Ramos and Taieb Touati (2015).
#' sensitivity: Sensitivity Analysis. R package version 1.11.1.
#' https://CRAN.R-project.org/package=sensitivity
#' 
#' @param n The number of samples
#' @param factors The model's parameters which will be evaluated
#' @param o Maximum order in the ANOVA decomposition
#' @param nb Number of bootstrap replicates
#' @param fun.doe The sampling function to be used for sobol method
#' @param fun.sobol The sobol implementation
#' 
#' 
#' @importFrom sensitivity sobol sobolmartinez
#' @export
AoE.Sobol<- function(n=100, factors=c(), o=2, nb=100, fun.doe=AoE.LatinHypercube, fun.sobol=sobolmartinez) {
  p.x1<- fun.doe(n,factors)
  p.x2<- fun.doe(n,factors)
  v<- fun.sobol(model = NULL, X1 = p.x1,X2 = p.x2, order = o, nboot = nb) 
  return(v)
}
  
##
## ----- Below Plot functions
##

#' @title Plot stability of output
#' 
#' @description Generate plot for visually access the stability of 
#' coefficient of variation as function of simulation sample size.
#' 
#' @param obj An instance of Morris Object \code{\link{AoE.Morris}}
#' @param title Chart title, may be null
#' 
#' @return The resulting ggplot2 plot object
#' 
#' @importFrom ggplot2 ggplot aes geom_point geom_line ggtitle labs geom_bar geom_errorbar
#' @export
Plot.Stability<- function(obj, title= NULL) {
  
  if(is.null(obj$RSD)) {
    stop("Invalid object instance!")
  }
  
  d<- obj
  
  p<- ggplot(d, with(d,aes(sample, RSD))) 
  
  p<- p + labs(y = expression("RSD"))
  p<- p + labs(x = expression("sample size"))
  
  if(!is.null(title)) {
    p<- p + ggtitle(title)
  }
  
  p<- p + with(d,aes(shape = group)) 
  p<- p + geom_line( with(d,aes(colour = group)), size = 1)
  
  return(p)
}

#' @title Plot of Morris output
#' 
#' @description Generate plot for Morris's screening method
#' 
#' @param obj An instance of Morris Object \code{\link{AoE.Morris}}
#' @param type The chart type (mu*sigma|musigma|mu*mu)
#' @param title Chart title, may be null
#' 
#' @return The resulting ggplot2 plot object
#' 
#' @importFrom ggplot2 ggplot aes geom_point geom_line ggtitle labs geom_bar geom_errorbar
#' @export
Plot.Morris<- function(obj, type, title= NULL) {
  # --- Check if we received a valid morris object
  if(is.null(obj$call)) {
    stop("Invalid Morris object instance!")
  }
  
  d<- AoE.GetMorrisOutput(obj)
  
  switch(type,
    "mu*sigma" = { 
      p<- ggplot(d, with(d,aes(mu.star, sigma)))  
      p<- p + labs(y = expression(sigma))
      p<- p + labs(x = expression(paste(mu,"*")))
    },
    
    "musigma" = { 
      p<- ggplot(d, with(d,aes(mu, sigma)))
      p<- p + labs(y = expression(sigma))
      p<- p + labs(x = expression(mu))
    },
    
    "mu*mu" = {
      p<- ggplot(d, with(d,aes(mu.star, mu)))
      p<- p + labs(y = expression(mu))
      p<- p + labs(x = expression(paste(mu,"*")))
    },
    
    stop("Invalid chart type!")
  )
  
  if(!is.null(title)) {
    p<- p + ggtitle(title)
  }
  
  p<- p + with(d,aes(shape = group)) 
  p<- p + geom_point( with(d,aes(colour = group)), size = 4)
  p<- p + geom_point(colour="grey90", size = 1.5)
  
  return(p)
}

#' @title Plot of Sobol output
#' 
#' @description Generate plot for Sobol's GSA
#' 
#' @param obj An instance of Sobol Object \code{\link{AoE.Sobol}}
#' @param type The chart type 
#' @param title Chart title, may be null
#' 
#' @return The resulting ggplot2 plot object
#' 
#' @importFrom ggplot2 ggplot aes geom_point geom_line ggtitle labs geom_bar geom_errorbar
#' @export
Plot.Sobol<- function(obj, type, title= NULL) {
  # --- Check if we received a valid sobol object
  if(is.null(obj$S)) {
    stop("Invalid Sobol object instance!")
  }
  
  switch(type,
         # --- First order indices
         "1" = { 
           d<- obj$S
           # --- Add the group column based on rownames
           d$group<- rownames(obj$S)
           
           y.label<- labs(y = expression(S[i]))
         },
         # --- Total order indices
         "2" = { 
           d<- obj$T
           # --- Add the group column based on rownames
           d$group<- rownames(obj$T)
           
           y.label<- labs(y = expression(S[Ti]))
         },
         
         stop("Invalid chart type!")
  )
  
  # --- Create plot object
  p<- ggplot(d, with(d,aes(group,original)))
  p<- p + y.label
  p<- p + labs(x = expression(paste("parameter")))

  ## --- p<- p + geom_bar(stat="identity",aes(fill = group))
  p<- p + geom_bar(stat="identity")
  p<- p + geom_errorbar( with(d, aes(ymin=`min. c.i.`, ymax=`max. c.i.`)), colour="black", width=.1)
  
  if(!is.null(title)) {
    p<- p + ggtitle(title)
  }
  
  return(p)
}
 
#' @title Plot of calibration 
#' 
#' @description Generate plot for parameter sets 
#' providing best fit
#' 
#' @param obj An instance of calibration Object
#' @param key The column name
#' @param title Chart title, may be null
#' 
#' @return The resulting ggplot2 plot object
#' 
#' @importFrom ggplot2 ggplot aes geom_point geom_line ggtitle labs geom_bar geom_errorbar
#' @export
Plot.Calibration<- function(obj, key, title= NULL) {
  # --- Check if we received a valid sobol object
  if(is.null(obj$pset)) {
    stop("Invalid calibration object instance")
  }
  
  d<- obj
  d$pset<- factor(d$pset, levels = d$pset)
  
  y.label<- labs(y = expression("Goodness of fit"))

  # --- Create plot object
  p<- ggplot(d, with(d,aes(pset,d[,key])))
  p<- p + y.label
  p<- p + labs(x = expression(paste("Parameter set")))
  
  ## --- p<- p + geom_bar(stat="identity",aes(fill = group))
  p<- p + geom_bar(stat="identity")

  if(!is.null(title)) {
    p<- p + ggtitle(title)
  }
  
  return(p)
}

##
## ----- Below Easy Api Methods
##

#' @title Easy.getChart
#' 
#' @description Returns the chart instance
#' 
#' @param obj A reference to the output of Easy.Stability
#' @param key The param name
#' 
#' @return The plot instance
#' @export
Easy.getChart<- function(obj, key) {
  if(is.null(obj$charts)) {
    stop("Not an instance of Easy API result!")
  }
  charts<- obj$charts
  chart<- charts[charts[,1] ==  key,]
  return(chart)
}

#' @title Easy API for output stability
#' 
#' @description This functions run model several times in order to determine 
#' how many experiment replications are required for model's output being stable
#' (i.e. the convergence of standard deviation)
#' 
#' @param m.dir The installation directory of some repast model
#' @param m.ds The name of any model aggregate dataset
#' @param m.time The total simulated time
#' @param parameters The factors or model's parameter list
#' @param samples The number of factor samples.
#' @param tries The number of experiment replications
#' @param vars The model's output variables for compute CoV
#' @param FUN The calibration function.
#' 
#' @return A list with holding experimnt, object and charts 
#' 
#' @export
Easy.Stability<- function(m.dir, m.ds, m.time=300, parameters, samples=1, tries=100, vars= c(), FUN) {
  my.model<- Model(modeldir=m.dir,maxtime = m.time, dataset=m.ds)
  Load(my.model)

  ## --- Sample the parameter space
  sampling<- AoE.RandomSampling(samples, parameters)
  
  ## --- Get the model declared paramters
  parms<- GetSimulationParameters(my.model)
  
  ## --- Build the experimental parameter set
  exp.design<- BuildParameterSet(sampling, parms)
  
  ## --- Run the experimental setup
  exp<- RunExperiment(my.model,r=tries,exp.design,FUN)
  
  ## --- Get the raw data set for evaluate the Coefficient of Variation
  d<- getExperimentDataset(exp)
  
  ## --- Calculate the coefficient of variation
  rsd<- AoE.Stability(d, vars)
  
  charts<- c()
  for(group in unique(rsd$group)) {
    chart<- Plot.Stability(rsd[rsd$group == group, ],"Simulation output stability")  
    charts<- rbind(charts, list(group=group,plot=chart))
  }
  
  if(length(vars) != 0) {
    chart<- Plot.Stability(rsd,"Simulation output stability")
    charts<- rbind(charts, list(group="all",plot=chart))
  }
  
  results<- list(experiment=exp, object=rsd, charts=charts)
  return(results)
  
}

#' @title Easy API for Morris's screening method
#' 
#' @description This functions wraps all calls to perform Morris method.
#' 
#' @param m.dir The installation directory of some repast model
#' @param m.ds The name of any model aggregate dataset
#' @param m.time The total simulated time
#' @param parameters The factors for morris screening.
#' @param mo.p The number of levels for the model's factors.
#' @param mo.r Repetitions. The number of random sampling points of Morris Method.
#' @param exp.r The number of experiment replications
#' @param FUN The calibration function.
#' 
#' @return A list with holding experimnt, object and charts 
#' 
#' @importFrom sensitivity tell
#' 
#' @export
Easy.Morris<- function(m.dir, m.ds, m.time=300, parameters, mo.p, mo.r, exp.r, FUN) {
  my.model<- Model(modeldir=m.dir,maxtime = m.time, dataset=m.ds)
  Load(my.model)
  
  ## --- Create Morris object
  v.morris<- AoE.Morris(parameters,p=mo.p,r=mo.r)
  
  ## --- Get the model declared paramters
  parms<- GetSimulationParameters(my.model)
  
  ## --- Build the experimental parameter set
  exp.design<- BuildParameterSet(v.morris$X,parms)
  
  ## --- Run the experimental setup
  exp<- RunExperiment(my.model,r=exp.r,exp.design,FUN)
  
  charts<- c()
  o<- getExperimentOutput(exp)
  for(k in colnames(o)) {
    if(k != "pset") {
      m<- t(df2matrix(getExperimentOutput(exp),c(k)))
      tell(v.morris,m)
      
      ## --- Plot Morris output
      mustar<- Plot.Morris(v.morris,"mu*sigma", paste("criteria",k))
      musigma<- Plot.Morris(v.morris,"musigma", paste("criteria",k))
      mumu<- Plot.Morris(v.morris,"mu*mu", paste("criteria",k))
      charts<- rbind(charts,list(mu.star=mustar,mu=musigma,mumu=mumu))
    } 
    ### ---> results<- list(experiment=exp, object=v.morris, charts=charts)
  }
  
  results<- list(experiment=exp, object=v.morris, charts=charts)
  return(results)
}

#' @title Easy API for Sobol's SA method
#' 
#' @description This functions wraps all required calls to perform 
#' Sobol method for global sensitivity analysis.
#' 
#' @param m.dir The installation directory of some repast model
#' @param m.ds The name of any model aggregate dataset
#' @param m.time The total simulated time
#' @param parameters The input factors
#' @param exp.n The experiment sample size
#' @param exp.r The number of experiment replications
#' @param bs.size The bootstrap sample size for sobol method
#' @param FUN The calibration function.
#' 
#' @return A list with holding experimnt, object and charts 
#' 
#' @importFrom sensitivity tell
#' 
#' @export
Easy.Sobol<- function(m.dir, m.ds, m.time=300, parameters,exp.n = 500, bs.size = 200, exp.r=1, FUN) {
  ## --- Instantiate the model
  my.model<- Model(modeldir=m.dir,maxtime = m.time, dataset=m.ds)
  Load(my.model)
  
  ## --- Get the model declared paramters
  parms<- GetSimulationParameters(my.model)
  
  ## --- Create a Sobol object
  my.obj<- AoE.Sobol(n= exp.n, parameters, nb=bs.size)
  
  # Build the experimental parameter set
  exp.design<- BuildParameterSet(my.obj$X,parms)
  
  ## --- Run the experimental setup
  exp<- RunExperiment(my.model,r=exp.r,exp.design,FUN)
  
  charts<- c()
  o<- getExperimentOutput(exp)
  for(k in colnames(o)) {
    if(k != "pset") {
      m<- t(df2matrix(getExperimentOutput(exp),c(k)))
      tell(my.obj,m)
      
      # -- First order indexes
      chart_0<- Plot.Sobol(my.obj, 1, paste("Sobol indexes for", k))
      
      # -- Total order indexes
      chart_1<- Plot.Sobol(my.obj, 2, paste("Sobol indexes for", k))
      
      charts<- rbind(charts,list(chart=chart_0))
      charts<- rbind(charts,list(chart=chart_1))
    } 
    ### ---> results<- list(experiment=exp, object=my.obj, charts=charts)
  }
  results<- list(experiment=exp, object=my.obj, charts=charts)
  return(results)
}

#' @title Easy.Setup
#' 
#' @description This function configures the deployment directory 
#' where logs and output dataset will be generated.  By default 
#' the deployment directory will be created under the model 
#' installation directory. The output generated by the Repast model 
#' will be redirected to the SystemOut.log file.  
#' 
#' @details If the deployment directory is empty the installation 
#' directory given by the parameter \code{model} is used instead as 
#' the base directory. The deployment directory is \code{/rrepast-deployment/}.
#' 
#' @param model The base directory where Repast model is installed.
#' @param deployment The directory to save the output and logs.
#' 
#' @export
Easy.Setup<- function(model, deployment=c()){
  
  if(length(deployment) == 0) {
    deployment<- paste0(model,"/rrepast-deployment/")
  }
  
  setOutputDir(deployment)
  
  ## -- Create output dir if required
  createOutputDir()
  
  jvm.setOut("SystemOut.log")
  PB.enable()
}

#' @title Easy.Calibration
#' 
#' @description Search for the best set of parameters trying to 
#' minimize the calibration function provided by the user. The function 
#' has to operational models, the first based on the experimental setup 
#' where all parameters are defined a priori and the second using 
#' optimization techniques. Currently the only supported optimization 
#' technique is the particle swarm optimization.
#' 
#' @param m.dir The installation directory of some repast model
#' @param m.ds The name of any model aggregate dataset
#' @param m.time The total simulated time
#' @param parameters The input factors
#' @param exp.n The experiment sample size
#' @param exp.r The number of experiment replications
#' @param smax The number of solutions to be generated
#' @param design The sampling scheme ["lhs"|"mcs"|"ffs"]
#' @param FUN The calibration function.
#'
#' @return A list with holding experiment, object and charts 
#' 
#' @examples \dontrun{
#'  my.cost<- function(params, results) {
#'    criteria<- c()
#'    Rate<- AoE.RMSD(results$X.Simulated,results$X.Experimental)
#'    G<- AoE.RMSD(results$G.T.,52)
#'    total<- Rate + G
#'    criteria<- cbind(total,Rate,G)
#'    return(criteria)
#'  }
#'  
#'  Easy.Setup("/models/BactoSim")
#'  v<- Easy.Calibration("/models/BactoSim","ds::Output",360,
#'                        f,exp.n = 1000, exp.r=1, smax=4, 
#'                        design="mcs", my.cost)
#'  
#' }
#' 
#' @export
Easy.Calibration<- function(m.dir, m.ds, m.time=300, parameters, exp.n = 100, exp.r=1, smax=4, design="lhs", FUN) {
  ## --- Sample the parameter space
  
  method<- "simple"
  switch(method,
    simple = {
      v<- simple.fitting(m.dir, m.ds, m.time, parameters, exp.n, exp.r, design , smax, FUN)      
    }, 
    stop("Valid calibration methods are [simple|pso]")
  )
  
  return(v)
}


##
## ----- Below calibration support methods
##


#' @title simple.fitting
#' 
#' @description Simple calibration method. Run an experimental setup and select the 
#' the best results minimizing the calibration function
#' 
#' @param m.dir The installation directory of some repast model
#' @param m.ds The name of any model aggregate dataset
#' @param m.time The total simulated time
#' @param parameters The input factors
#' @param samples The experiment sample size
#' @param tries The number of experiment replications
#' @param design The sampling scheme ["lhs"|"mcs"|"ffs"]
#' @param smax The number of solutions to be generated
#' @param objective The calibration function.
#'
#' @importFrom gridExtra ttheme_default tableGrob arrangeGrob
#' @export
simple.fitting<- function(m.dir, m.ds, m.time=300, parameters, samples=100, tries=1, design="lhs" , smax=4, objective) {
  ## --- Instantiate the model
  my.model<- Model(modeldir=m.dir,maxtime = m.time, dataset=m.ds,load = TRUE)
  
  
  ## --- Sample the parameter space
  switch(design,
    lhs = {
      sampling<- AoE.LatinHypercube(samples, parameters)  
    },
    
    mcs = {
      sampling<- AoE.RandomSampling(samples, parameters)
    },
    
    ffs = {
      sampling<- AoE.FullFactorial(samples, parameters)
    },
    
    stop("Valid sampling types are [mcs|lhs|ffs]")
    
  )
         
  
  ## --- Get the model declared paramters
  parms<- GetSimulationParameters(my.model)
  
  ## --- Build the experimental parameter set
  exp.design<- BuildParameterSet(sampling, parms)
  
  ## --- Run the experimental setup
  exp<- RunExperiment(my.model,r=tries,exp.design,objective)
  
  ## --- Add a totalization column
  exp$output<- col.sum(exp$output)

  tbl.theme<- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
  
  ##tmp<- c()
  charts<- c()
  obj<- c()
  fittest.max<- smax 
  
  o<- getExperimentOutput(exp)
  for(k in colnames(o)) {
    if(k != "pset") {
      p<- getExperimentParamSet(exp)
      best<- dffilterby(o,"pset",pick.fittest(o,goals=c(k),fittest.max)$pset)
      chart<- Plot.Calibration(best,k,paste0("Best parameters for ", k))
      
      best.p<- dffilterby(p,"pset",pick.fittest(o,goals=c(k),fittest.max)$pset)
      tbl.data<- best.p[,c("pset",parameters[,"name"])]
      tbl.data<- dfround(tbl.data,2)
      tbl.table<- tableGrob(tbl.data, rows=NULL, theme= tbl.theme)
      my.chart<- arrangeGrob(chart, tbl.table, nrow=2, as.table=TRUE, heights=c(3,1))
      
      charts<- rbind(charts,list(variable=k,both=my.chart,chart=chart,table=tbl.table))
      obj<- rbind(obj,list(variable=k,parameters=best.p,objective=best))
    }
  }
  
  ###obj$data<- tmp
  ###obj$keys<- unlist(obj$data[,"variable"])
  ###obj$parameters<- function(k) {obj$data[(obj$data[,"variable"] == k),"parameters"][[1]]}
  ###obj$objective<- function(k) {obj$data[(obj$data[,"variable"] == k),"objective"][[1]]}

  results<- list(experiment=exp, object=obj, charts=charts)
  return(results)
}

##
## ----- Functions for accessing result object members
##

#' @title Results.GetExperiment
#' 
#' @description Simplify the access to the experiment member
#' 
#' @param obj An instance of the object returned by \code{Easy} methods
#' 
#' @return The experiment element inside results
#' @export
Results.GetExperiment<- function(obj) {
  if(is.null(obj$experiment)) {
    stop("Not an instance of Easy API result!")
  }
  obj$experiment
}

#' @title Results.GetObject
#' 
#' @description Simplify the access to the object member
#' 
#' @param obj An instance of the object returned by \code{Easy} methods
#' 
#' @return The object element inside results
#' @export
Results.GetObject<- function(obj) {
  if(is.null(obj$object)) {
    stop("Not an instance of Easy API result!")
  }
  obj$object
}

#' @title Results.GetCharts
#' 
#' @description Simplify the access to the charts member
#' 
#' @param obj An instance of the object returned by \code{Easy} methods
#' 
#' @return The charts element inside results
#' @export
Results.GetCharts<- function(obj) {
  if(is.null(obj$charts)) {
    stop("Not an instance of Easy API result!")
  }
  obj$charts
}

#' @title Calibration.GetMemberKeys
#' 
#' @description Gets the list of keys (the factor names)
#' 
#' @param obj An instance of the object returned by \code{Easy} methods
#' 
#' @return The collection of keys
#' @export
Calibration.GetMemberKeys<- function(obj) {
  if(!"variable" %in% colnames(obj)) {
    stop("Not an instance of a Easy.Calibration return!")
  }
  unlist(obj[,"variable"])
}

#' @title Calibration.GetMemberList
#' 
#' @description Gets the member list value
#' 
#' @param obj An instance of the object returned by \code{Easy} methods
#' @param key The key value
#' @param name The column name 
#' 
#' @return The member list
#' @export
Calibration.GetMemberList<- function(obj, key, name) {
  if(!"variable" %in% colnames(obj)) {
    stop("Not an instance of a Easy.Calibration return!")
  }
  obj[(obj[,"variable"] == key),name][[1]]
}


#' @title pick.fittest
#' 
#' @description Choose the best solutions minimizing the objective function
#' 
#' @param out The output data set holding the values of goals 
#' @param goals The column names which must be used as goal
#' @param n The number of solutions
#' 
#' @return The n rows holding the best results
#' 
#' @export
pick.fittest<- function(out, goals=c(), n=4) {
  out<- as.data.frame(out) 
  
  ## -- Check if out was generated by RunExperiment
  if(!"pset" %in% colnames(out)) {
    stop("Invalid data set!")  
  }
  
  ## --- Adjusting defaults
  n<- ifelse(n > nrow(out),nrow(out),n)
  goals<- ifelse(length(goals) == 0,c(2),goals)
  
  out<- out[order(out[,goals]), ]
  return(out[1:n,])
}

#' @title col.sum
#' 
#' @description Sum all columns but one (pset) of a data frame
#' 
#' @param d The data frame 
#' @param skip The columns which should not be included in the sum
#' 
#' @return The original data frame with a new column (sum) holding the sum 
#' 
#' @export
col.sum<- function(d,skip=c()) {
  v<- as.data.frame(d)
  s<- NULL
  
  ## -- always skip pset
  if(!"pset" %in% skip) {
    skip<- c("pset",skip)  
  }
  
  ## -- always skip total
  if(!"total" %in% skip) {
    skip<- c("total",skip)  
  }
    
  v<- dfsumcol(v,skip)
  
  return(v)
}

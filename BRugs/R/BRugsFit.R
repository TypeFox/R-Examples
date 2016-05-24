BRugsFit <-
function(modelFile, data, inits, numChains = 3, parametersToSave,
    nBurnin = 1000, nIter = 1000, nThin = 1, coda = FALSE,
    DIC = TRUE, working.directory = NULL, digits = 5, seed=NULL,
    BRugsVerbose = getOption("BRugsVerbose")){

  if(is.null(BRugsVerbose))
      BRugsVerbose <- TRUE
  op <- options("BRugsVerbose" = BRugsVerbose)
  on.exit(options(op))
  if(!is.null(working.directory)){
      working.directory <- path.expand(working.directory)
      savedWD <- getwd()
      setwd(working.directory)
      on.exit(setwd(savedWD), add = TRUE)
  }
  if(is.function(modelFile)){
      writeModel(modelFile, con = (modelFile <- tempfile("model")), digits = digits)
      if(!is.R()) on.exit(file.remove(modelFile), add = TRUE)
  }
  if(!file.exists(modelFile)) stop(modelFile, " does not exist")
  if(file.info(modelFile)$isdir) stop(modelFile, " is a directory, but a file is required")
  modelCheck(modelFile)
  if(!(is.vector(data) && is.character(data) && all(file.exists(data))))
    data <- bugsData(data, digits = digits)
  modelData(data)
  modelCompile(numChains)
  if(!is.null(seed)) modelSetRN(seed)
  if(!missing(inits)){
      if(is.list(inits) || is.function(inits))
          inits <- bugsInits(inits = inits, numChains = numChains, digits = digits)
      if (is.character(inits) && any(file.exists(inits))){
          if(BRugsVerbose) print(inits)
          modelInits(inits)
      }
  }
  modelGenInits()
  samplesSetThin(nThin)
  modelUpdate(nBurnin)
  if(DIC){
    dicSet()
    on.exit(dicClear(), add = TRUE)
  }
  samplesSet(parametersToSave)
  modelUpdate(nIter)
  if(coda)
    return(buildMCMC("*"))
  else
    return(list(Stats = samplesStats("*"), DIC = if(DIC) dicStats()))
}

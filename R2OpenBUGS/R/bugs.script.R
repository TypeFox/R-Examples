"bugs.script" <-
  function(parameters.to.save, n.chains, n.iter, n.burnin,
           n.thin, saveExec, restart, model.file.bug,
           model.file, debug=FALSE, is.inits, 
           DIC=FALSE, useWINE=FALSE,
           newWINE=TRUE, WINEPATH=NULL, bugs.seed=NULL, summary.only=FALSE,
           save.history=(.Platform$OS.type == "windows" | useWINE==TRUE),
           bugs.data.file, bugs.inits.files,
           over.relax = FALSE)
{
  ## Write file script.txt for Bugs
  if(n.iter - n.burnin < 2)
    stop ("(n.iter-n.burnin) must be at least 2")
  working.directory <- getwd()
  script <- "script.txt"

  model <- 
    if (length(grep("\\\\", model.file)) || length(grep("/", model.file))) {
        gsub("\\\\", "/", model.file)
    }
    else file.path(working.directory, model.file)
  model <- native2win(model, useWINE=useWINE, newWINE=newWINE, WINEPATH=WINEPATH)

  data <- file.path(working.directory, bugs.data.file)
  data <- native2win(data, useWINE=useWINE, newWINE=newWINE, WINEPATH=WINEPATH)

  coda  <- file.path(working.directory, "/")
  coda <- native2win(coda, useWINE=useWINE, newWINE=newWINE, WINEPATH=WINEPATH)

  model.file.bug<-file.path(working.directory,model.file.bug)
  model.file.bug<-native2win(model.file.bug, useWINE=useWINE, newWINE=newWINE, WINEPATH=WINEPATH)

  logFile <- file.path(working.directory, "log.odc")
  logFile <- native2win(logFile, useWINE=useWINE, newWINE=newWINE, WINEPATH=WINEPATH)
  logFileTxt <- file.path(working.directory, "log.txt")
  logFileTxt <- native2win(logFileTxt, useWINE=useWINE, newWINE=newWINE, WINEPATH=WINEPATH)

  inits <- paste(working.directory, "/", bugs.inits.files, sep="")
  inits <- sapply(inits, useWINE=useWINE, newWINE=newWINE, WINEPATH=WINEPATH, 
    function(x, useWINE, newWINE, WINEPATH) 
    {native2win(x, useWINE=useWINE, newWINE=newWINE, WINEPATH=WINEPATH)})

  initlist <- paste("modelInits(", "'", inits, "',",1:n.chains,")\n", sep="")

  savelist <- paste("samplesSet(", parameters.to.save, ")\n", sep="")
  summarylist <- paste("summarySet(", parameters.to.save, ")\n", sep="")

  bugs.seed.cmd <- ""
  if (!is.null(bugs.seed)) {
        bugs.seed.cmd <- paste("modelSetRN(", bugs.seed, ")\n", sep="")
  }
    
  thinUpdate <- paste("modelUpdate(", formatC(n.burnin, format='d'), ",", n.thin, 
                      ",",formatC(n.burnin, format='d'), ")\n", sep="")

  cat(
    if(.Platform$OS.type == "windows" | useWINE) "modelDisplay('log')\n",
    if(restart)c("modelInternalize('", model.file.bug, "')\n"),
    if(restart && n.burnin>0)c(
       "samplesClear('*')\n",
       "summaryClear('*')\n"
       ),
    if(!restart)c( 
      "modelCheck('", model, "')\n",
      "modelData('", data, "')\n",
      "modelCompile(", n.chains, ")\n"
      ),
    if(!restart)bugs.seed.cmd,
    if(!restart && is.inits) initlist,
    if(!restart)"modelGenInits()\n",
    if(!restart && over.relax) 'over.relax("yes")\n',
    if((!restart) || (n.burnin>0))c(
    thinUpdate,
    savelist,
    summarylist
    ),
    if(((!restart) || (n.burnin>0)) && DIC) "dicSet()\n",
    "modelUpdate(", formatC(n.iter-n.burnin, format='d'), ",", n.thin, 
                  ",",formatC(n.iter-n.burnin, format='d'),")\n", 
    "samplesCoda('*', '", coda, "')\n", 
    "summaryStats('*')\n",
    if(DIC) "dicStats()\n",
    if (save.history) "samplesHistory('*')\n", 
    if(saveExec)c("modelExternalize('",model.file.bug,"')\n"),
    if(.Platform$OS.type == "windows" | useWINE) c("modelSaveLog('", logFile, "')\n",
    "modelSaveLog('", logFileTxt, "')\n"),
    file=script, sep="", append=FALSE)

    if(!debug) cat("modelQuit('y')\n", file=script, append=TRUE)

  sims.files <- paste("CODAchain", 1:n.chains, ".txt", sep="")
  for(i in 1:n.chains)
    cat("OpenBUGS did not run correctly.\n", file=sims.files[i], append=FALSE)
}

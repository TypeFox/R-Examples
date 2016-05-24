"bugs.script" <-
  function(parameters.to.save, n.chains, n.iter, n.burnin,
           n.thin, model.file, debug=FALSE, is.inits, bin,
           DIC=FALSE, useWINE=.Platform$OS.type != "windows",
           newWINE=TRUE, WINEPATH=NULL, bugs.seed=NULL, summary.only=FALSE,
           save.history=TRUE, bugs.data.file, bugs.inits.files,
           over.relax = FALSE)
{
  ## Write file script.txt for Bugs
  if((ceiling(n.iter/n.thin) - ceiling(n.burnin/n.thin)) < 2)
    stop ("(n.iter-n.burnin)/n.thin must be at least 2")
  working.directory <- getwd()
  script <- "script.txt"

  model <- 
    if (length(grep("\\\\", model.file)) || length(grep("/", model.file))) {
        gsub("\\\\", "/", model.file)
    }
    else file.path(working.directory, model.file)
  model <- native2win(model, useWINE=useWINE, newWINE=newWINE, WINEPATH=WINEPATH)

  history <- file.path(working.directory, "history.odc")
  history <- native2win(history, useWINE=useWINE, newWINE=newWINE, WINEPATH=WINEPATH)
  
  data <- file.path(working.directory, bugs.data.file)
  data <- native2win(data, useWINE=useWINE, newWINE=newWINE, WINEPATH=WINEPATH)

  coda  <- file.path(working.directory, "coda")
  coda <- native2win(coda, useWINE=useWINE, newWINE=newWINE, WINEPATH=WINEPATH)

  logFile <- file.path(working.directory, "log.odc")
  logFile <- native2win(logFile, useWINE=useWINE, newWINE=newWINE, WINEPATH=WINEPATH)
  logFileTxt <- file.path(working.directory, "log.txt")
  logFileTxt <- native2win(logFileTxt, useWINE=useWINE, newWINE=newWINE, WINEPATH=WINEPATH)

  inits <- paste(working.directory, "/", bugs.inits.files, sep="")
  inits <- sapply(inits, useWINE=useWINE, newWINE=newWINE, WINEPATH=WINEPATH, 
    function(x, useWINE, newWINE, WINEPATH) 
    {native2win(x, useWINE=useWINE, newWINE=newWINE, WINEPATH=WINEPATH)})

  initlist <- paste("inits (", 1:n.chains, ", '", inits, "')\n", sep="")

  savelist <- paste("set (", parameters.to.save, ")\n", sep="")
  redo <- ceiling((n.iter-n.burnin)/(n.thin*bin))

  bugs.seed.cmd <- ""
  if (!is.null(bugs.seed)) {
        bugs.seed.cmd <- paste("set.seed(", bugs.seed, ")\n", sep="")
  }
    
  thinUpdate <- paste("thin.updater (", n.thin, ")\n",
                      "update (", ceiling(n.burnin/n.thin), ")\n", sep="")

  cat(
    "display ('log')\n",
    "check ('", model, "')\n",
    "data ('", data, "')\n",
    bugs.seed.cmd,
    "compile (", n.chains, ")\n",
    if(is.inits) initlist,
    "gen.inits()\n",
    if(over.relax) 'over.relax("yes")\n',
    thinUpdate,
    savelist,
    if(DIC) "dic.set()\n",
    rep(c("update (", formatC(ceiling(bin), format = "d"), ")\n",
          #if (!summary.only) ## Hmm, if coda files are not written, we do not know if WinBUGS did not fail
          c("coda (*, '", coda, "')\n")), 
        redo),
    "stats (*)\n",
    if(DIC) "dic.stats()\n",
    if (save.history) c("history (*, '", history, "')\n"), 
    "save ('", logFile, "')\n",
    "save ('", logFileTxt, "')\n",
    file=script, sep="", append=FALSE)

  if(!debug) cat("quit ()\n", file=script, append=TRUE)

  sims.files <- paste("coda", 1:n.chains, ".txt", sep="")
  for(i in 1:n.chains)
    cat("WinBUGS did not run correctly.\n", file=sims.files[i], append=FALSE)
}

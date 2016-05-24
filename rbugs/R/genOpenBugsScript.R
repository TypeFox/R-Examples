#### OpenBugs Script Commands
ScriptCommands <- function(hateWindows = TRUE) {
  commands <- c("CHECK", "DATA", "COMPILE", "INITS",
                "GENINITS", "BEG", "UPDATE",
                "SET", "DICSET",
                "STATS", "DICSTATS", "CODA", "SAVE",
                "SETRN", "GETRN",
                "QUIT", "LBR")
  openBugs <- c("modelCheck", "modelData","modelCompile","modelInits",
                "modelGenInits", "samplesBeg", "modelUpdate",
                "samplesSet","dicSet",
                "samplesStats", "dicStats", "samplesCoda", "modelSaveLog",
                "modelSetRN", "modelGetRN",
                "modelQuit", "\n")
  winBugs <- c("check", "data", "compile", "inits",
               "gen.inits", "beg", "update",
               "set", "dic.set",
               "stats", "dic.stats", "coda", "save",
               "set.seed", "get.seed",
               "quit", "\n")
  comm <- if(hateWindows) openBugs else winBugs
  names(comm) <- commands
  comm
}

genBugsScript <-
  function(paramSet,
           n.chains,
           n.iter,
           n.burnin,
           n.thin,
           dic = FALSE,
           model.file,
           data.file,
           inits.files,
           bugsWorkingDir=getwd(), ## needs to be readable for BUGS
           script, #output
           debug=FALSE,
           OpenBugs=TRUE,
           Windows=TRUE, ## Modified by Marcos
           seed=NULL ## This number cannot be < -148 or >148. How strange!
           ) {
  if (n.chains != length(inits.files)) stop("length(inits.files) should equal n.chains.")
  ## n.iter <- n.burnin + n.thin * n.keep

  if(length(seed)==0) ## Modified by Marcos
        seed = floor(runif(1,-148,149)) ## Modified by Marcos
 
  ## add deviance to the paramSet list
  paramSet <- c(paramSet, "deviance")

  ## setup workingDir
  workingDir <- bugsWorkingDir

  #if (is.null(workingDir)) {
    #if (useWine) workingDir <- driveTr(bugsWorkingDir, .DriveTable)
    #else workingDir <- bugsWorkingDir
  #}

  #if (OpenBugs) useWine <- FALSE
  ## necessary if useWine == TRUE
  #if (useWine) {
  #  model.file <- sub(workingDir, bugsWorkingDir, model.file)
  #  data.file <- sub(workingDir, bugsWorkingDir, data.file)
  #  for (i in 1:length(inits.files))
  #    inits.files[i] <- sub(workingDir, bugsWorkingDir, inits.files[i])
  #}

  ## attach the command list
  comm <- ScriptCommands(OpenBugs)
  LBR <- comm["LBR"]
  ## attach(comm)
  ## on.exit(detach(comm))
  
  ## setup some file names
  coda  <- file.path(bugsWorkingDir, "coda")
  ## logodc <- file.path(bugsWorkingDir, "log.odc")
  logfile <- file.path(bugsWorkingDir, "log.txt")
  ## note that the order or arguments to INITS are different
  ## in WinBUGS and OpenBUGS
  initlist <- if (OpenBugs) paste(comm["INITS"], "(", "'", inits.files, "', ", 1:n.chains, ")", LBR, sep="") else paste(comm["INITS"], "(", 1:n.chains, ", '", inits.files, "')", LBR, sep="")
  savelist <- paste(comm["SET"], "(", paramSet, ")", LBR, sep="")
  ## write out to script.txt
  nburn <- ceiling(n.burnin / n.thin)
  nsamp <- ceiling((n.iter - n.burnin) / n.thin)
  cat (
       if(Windows && OpenBugs) "modelDisplay ('log')\n",
       comm["CHECK"], "('", model.file, "')", LBR,
       comm["DATA"], "('", data.file, "')", LBR,
       comm["COMPILE"], "(", n.chains, ")", LBR,
       comm["SETRN"], "(", seed, ")", LBR,
       initlist,
       comm["GENINITS"], "()", LBR,
       comm["BEG"], "(", nburn + 1, ")", LBR,
       comm["UPDATE"], "(", nburn, ", ", n.thin, ")", LBR,
       savelist,
       if (dic) c(comm["DICSET"], "()", LBR),
       comm["UPDATE"], "(", nsamp, ", ", n.thin, ")", LBR,
       comm["STATS"], "('*')", LBR,
       if (dic) c(comm["DICSTATS"], "(*)", LBR),
       comm["CODA"], "('*', '", coda, "')", LBR,
       ## "save ('", logodc, "')\n", 
       ## comm["SAVE"], "('", logfile, "')", LBR,
       ## modelSaveLog is only available on windows.
       ##if (OpenBugs) c(comm["QUIT"], "()", LBR)
       if (!Windows) c(comm["QUIT"], "('yes')", LBR), ## Modified by Marcos
       ##else c("modelSaveLog", "('", logfile, "')", LBR),
       if (Windows) c(comm["SAVE"], "('", logfile, "')", LBR), ## Modified by Marcos
       ##if (Windows && OpenBugs && !debug) c(comm["QUIT"], "('yes')", LBR), ## Modified by Marcos
       if (Windows && OpenBugs) c(comm["QUIT"], "('yes')", LBR), ## Modified by Marcos
       file=script, sep="", append=FALSE)
  if (!debug && !OpenBugs) cat (comm["QUIT"], "()", LBR, sep="", file=script, append=TRUE) ## Modified by Marcos
}

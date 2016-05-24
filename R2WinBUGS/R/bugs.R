"bugs" <-
function(data, inits, parameters.to.save, model.file="model.bug",
    n.chains=3, n.iter=2000, n.burnin=floor(n.iter / 2),
    n.thin=max(1, floor(n.chains * (n.iter - n.burnin) / n.sims)), n.sims = 1000,
    bin=(n.iter - n.burnin) / n.thin,
    debug=FALSE, DIC=TRUE, digits=5, codaPkg=FALSE,
    bugs.directory="c:/Program Files/WinBUGS14/",
    program=c("WinBUGS", "OpenBUGS", "winbugs", "openbugs"),
    working.directory=NULL,
    clearWD=FALSE, useWINE=.Platform$OS.type != "windows", WINE=NULL,
    newWINE=TRUE, WINEPATH=NULL, bugs.seed=NULL, summary.only=FALSE,
    save.history=!summary.only, over.relax = FALSE)
{
  if(!is.null(working.directory)) {
    working.directory <- path.expand(working.directory)
    savedWD <- getwd()
    setwd(working.directory)
    on.exit(setwd(savedWD))
  }
  program <- match.arg(program)
  if (missing(bugs.directory) && 
        !is.null(bugs.dir <- getOption("R2WinBUGS.bugs.directory"))) { # requested by Jouni Kerman
      bugs.directory <- bugs.dir
  }
  if(program %in% c("openbugs", "OpenBUGS", "OpenBugs")) {
    if(!is.R()) stop("OpenBUGS is not yet available in S-PLUS")
    ## If OpenBUGS, we only call openbugs() and exit...
    return(openbugs(data, inits, parameters.to.save, model.file,
                    n.chains, n.iter, n.burnin, n.thin, n.sims, DIC=DIC,
                    bugs.directory, working.directory, digits, over.relax = over.relax, seed=bugs.seed))
  }
  ## Checking number of inits, which is NOT saved here:
  if(!missing(inits) && !is.function(inits) && !is.null(inits) && (length(inits) != n.chains))
    stop("Number of initialized chains (length(inits)) != n.chains")

  ## Wine
  if(useWINE) {
    if(!is.R())
      stop("Non-Windows platforms not yet supported in R2WinBUGS for S-PLUS")
    ## Attempt to find wine and winepath
    if(is.null(WINE)) WINE <- findUnixBinary(x="wine")
    if(is.null(WINEPATH)) WINEPATH <- findUnixBinary(x="winepath")
  }

  ## Move to working drirectory or temporary directory when NULL
  inTempDir <- FALSE
  if(is.null(working.directory)) {
    working.directory <- tempdir()
    if(useWINE){
        ## Some tweaks for wine (particularly required for Mac OS)
        working.directory <- gsub("//", "/", working.directory)
        Sys.chmod(working.directory, mode="770")
        on.exit(Sys.chmod(working.directory, mode="700"), add = TRUE)
    }
    savedWD <- getwd()
    setwd(working.directory)
    on.exit(setwd(savedWD), add = TRUE)
    inTempDir <- TRUE
  }

  ## model.file is not a file name but a model function
  if(is.function(model.file)){
      temp <- tempfile("model")
      temp <-
        if(is.R() || .Platform$OS.type != "windows"){
               paste(temp, "txt", sep=".")
        } else {
               gsub("\\.tmp$", ".txt", temp)
        }
      write.model(model.file, con=temp, digits=digits)
      model.file <- gsub("\\\\", "/", temp)
      if(!is.R()) on.exit(file.remove(model.file), add=TRUE)
  }
  if(inTempDir && basename(model.file) == model.file)
    try(file.copy(file.path(savedWD, model.file), model.file, overwrite = TRUE))
  if(!file.exists(model.file))
    stop(paste(model.file, "does not exist."))
  if(file.info(model.file)$isdir)
    stop(paste(model.file, "is a directory, but a file is required."))
  if (!(length(data) == 1 && is.vector(data) && is.character(data) && 
       (regexpr("\\.txt$", data) > 0))) {
    bugs.data.file <- bugs.data(data, dir = getwd(), digits)
  } else {
    if(inTempDir && all(basename(data) == data))
        try(file.copy(file.path(savedWD, data), data, overwrite = TRUE))
    if(!file.exists(data))
        stop("File", data, "does not exist.")
    bugs.data.file <- data
  }

  if (is.character(inits)) {
    if(inTempDir && all(basename(inits) == inits))
        try(file.copy(file.path(savedWD, inits), inits, overwrite = TRUE))
    if (!all(file.exists(inits))) {
        stop("One or more inits files are missing")
    }
    if (length(inits)!=n.chains) {
        stop("Need one inits file for each chain")
    }
    bugs.inits.files <- inits
  } else {
    if (!is.function(inits) && !is.null(inits) &&  (length(inits) != n.chains)) {
        stop("Number of initialized chains (length(inits)) != n.chains")
    }
    bugs.inits.files <- bugs.inits(inits, n.chains, digits)
  }

  if(DIC) parameters.to.save <- c(parameters.to.save, "deviance")
  ## Model files with extension ".bug" need to be renamed to ".txt"
  if(!length(grep("\\.txt$", tolower(model.file)))) {
    new.model.file <- paste(basename(model.file), ".txt", sep="")
    if(!is.null(working.directory)) new.model.file <- file.path(working.directory, new.model.file)
    file.copy(model.file, new.model.file, overwrite=TRUE)
    on.exit(try(file.remove(new.model.file)), add=TRUE)
  } else {
    new.model.file <- model.file
  }
  if(useWINE){
        ## Some tweaks for wine (particularly required for Mac OS)
        new.model.file <- gsub("//", "/", new.model.file)
  }
  bugs.script(parameters.to.save, n.chains, n.iter, n.burnin, n.thin,
              new.model.file, debug=debug, is.inits=!is.null(inits),
              bin=bin, DIC=DIC, useWINE=useWINE, newWINE=newWINE,
              WINEPATH=WINEPATH, bugs.seed=bugs.seed, 
              summary.only=summary.only, save.history=save.history, 
              bugs.data.file = bugs.data.file, 
              bugs.inits.files = bugs.inits.files, over.relax = over.relax)
  bugs.run(n.burnin, bugs.directory, WINE=WINE, useWINE=useWINE,
           newWINE=newWINE, WINEPATH=WINEPATH)
  if(codaPkg)
    return(file.path(getwd(), paste("coda", 1:n.chains, ".txt", sep="")))
  if (summary.only) {
    return(bugs.log("log.txt"))
  }

  sims <- c(bugs.sims(parameters.to.save, n.chains, n.iter, n.burnin,
                      n.thin, DIC),
            model.file=model.file, program=program)
  if(clearWD) {
    file.remove(c(bugs.data.file, "log.odc", "log.txt", "codaIndex.txt",
                  bugs.inits.files, "script.txt",
                  paste("coda", 1:n.chains, ".txt", sep="")))
  }
  class(sims) <- "bugs"
  sims
}

#### configurations  12/15/2003
## BUGS stores the executable of bugs
## workingDir is the directory to save all files, default is ??tempdir()??
## bugsWorkingDir is the directory for wine to use windows type directory
## WINE stores the executable of wine
## useWine = TRUE if use wine

## example usage on linux
## schools.sim <- rbugs(data=schools.data, inits, parameters, "schools.bug", n.chains=3, n.iter=1000, workingDir="/var/tmp/jyan/c/tmp", bugsWorkingDir="c:/tmp", useWine=T, wine="/var/scratch/jyan//wine/wine-20031016/wine", debug=T)

## Changed the linbugs

rbugs <- function(data, inits, paramSet, model,
                  ## mcmc options
                  n.chains=1, n.iter=2000, n.burnin=floor(n.iter/2),
                  n.thin=max(1, floor(n.chains*(n.iter-n.burnin)/1000)),
                  ## monitoring options
                  dic=FALSE,
                  ## configuration options
                  debug=FALSE,
                  bugs= system("which OpenBUGS", TRUE),
                  ##"c:/Program Files/WinBUGS14/WinBUGS14.exe",
                  #workingDir = NULL, #getwd(),
                  ##"/var/scratch/jyan/c/tmp", # native
                  bugsWorkingDir, # required argument
                  ##"c:/tmp",
                  OpenBugs = TRUE, ##Modified by Marcos ##Using OpenBugs or not???
                  cleanBugsWorkingDir = FALSE,
                  genFilesOnly = FALSE,
                  verbose = FALSE,
                  seed=NULL
                  ){
  ##  start.time <- Sys.time ()
  ##linbugs = OpenBugs ##Modified by Marcos
  Windows = FALSE ##Modified by Marcos
  os.type <- .Platform$OS.type
  if(length(bugsWorkingDir) == 0) stop("It is required to specify the bugsWorkingDir")
  if (os.type == "windows") {
    if (!file.exists(bugs))
      stop(paste("BUGS executable", bugs, "does not exists."))
    Windows = TRUE ##Modified by Marcos
  }
  else if (os.type == "unix") {
      if (length(bugs) == 0) bugs <- system("which OpenBUGS", TRUE)
      if (length(bugs) == 0)
        stop(paste("BUGS executable", bugs, "does not exists."))
      # bugs <- filePathAsAbsolute(bugs) ##Modified by Marcos 
  }
  else warning("This function has not been tested on mac-os.")
  
  ## setup workingDir
  #bugsWorkingDir <- filePathAsAbsolute(bugsWorkingDir)
  if (is.null(bugsWorkingDir)) {
    bugsWorkingDir <- getwd() #tempfile("bugsWorkingDir")
    if (!file.exists(bugsWorkingDir)) dir.create(bugsWorkingDir)
    on.exit(if(cleanBugsWorkingDir) unlink(bugsWorkingDir, TRUE))
  }
  workingDir <- bugsWorkingDir ## Modefified by Marcos
  #workingDir <- filePathAsAbsolute(workingDir)
  
  ## prepare the model file by 
  ## making a copy of model to the working directory
  if (!file.exists(model)) stop("Model file doesn't exits.")
  model.file <- file.path(workingDir, "model.txt")
  file.copy(model, model.file, overwrite=TRUE)

  ## prepare the data file
  data.file <- file.path(workingDir, "data.txt")
  genDataFile(data, data.file)

  ## prepare the inits files
  inits.file.stem <- file.path(workingDir, "init")
  genInitsFile(n.chains, inits, inits.file.stem)
  inits.files <- paste(inits.file.stem, 1:n.chains, ".txt", sep="")

  ## prepare the script file
  script.file <- paste(workingDir, "script.txt", sep="/")
  genBugsScript(paramSet, n.chains, n.iter, n.burnin, n.thin, dic,
                model.file, data.file, inits.files,
                bugsWorkingDir,
                script.file, debug, OpenBugs, Windows, seed) ##Modified by Marcos

  ## change line breaks from "\n" to "\r\n"
  ## otherwise, OpenBugs would hang!!
  if (OpenBugs) {
    ## trLbr(script.file)
    ## spend three hours to figure out that script file doesn't need it!
    trLbr(model.file)
    trLbr(data.file)
    for (i in inits.files) trLbr(i)
  }
  
  ## run bugs
  if (genFilesOnly) {
    cat("Files are generated in", workingDir, "\n")
    return(TRUE)
  }
  #if (useWine) script.file <- gsub(workingDir, bugsWorkingDir, script.file)
  #runBugs(bugs, script.file, n.chains, workingDir, useWine, wine, OpenBugs, Windows, verbose)
  runBugs(bugs, script.file, n.chains, workingDir, OpenBugs, Windows, verbose)

  ## collect the output
#  out <- getBugsOutput(n.chains, workingDir, OpenBugs)
#  all <- list (n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, ##Modified by Marcos
#               n.thin=n.thin, n.keep=out[[n.chains+1]], n.sims=out[[n.chains+2]]) ##Modified by Marcos
#  out <- out[-(n.chains+2)] ##Modified by Marcos
#  out <- out[-(n.chains+1)] ##Modified by Marcos
#  all <- c(all, out) ##Modified by Marcos 
#  class(all) <- "rbugs" ##Modified by Marcos

  all <- getBugsOutput(n.chains, workingDir, OpenBugs)  ##Modified by Marcos
  all$n.iter = n.iter                                   ##Modified by Marcos
  all$n.burnin = n.burnin                               ##Modified by Marcos
  all$n.thin = n.thin                                   ##Modified by Marcos
  
  if(cleanBugsWorkingDir) { ##Modified by Marcos
    fnames <- getCodaFileNames(n.chains, workingDir, OpenBugs) ##Modified by Marcos
    coda.files <- c(fnames$codaIndexFile, fnames$codaFiles) ##Modified by Marcos
    log.file <- file.path(workingDir, "log.txt") ##Modified by Marcos
    file.remove(c(model.file, log.file, data.file, ##Modified by Marcos
                  inits.files, script.file, coda.files)) ##Modified by Marcos
  }
 
  all
}


genDataFile <- function(dataList, dataFile) {
  if (is.numeric(unlist(dataList))) {
    ## cat(dput2bugs(dataList), file = data.file)
    cat(format4Bugs(dataList), file = dataFile, fill = TRUE)
  }
  else {
    data <- lapply(dataList, get, pos = 1)
    names(data) <- dataList
    ## cat(dput2bugs(data), file = data.file, fill = TRUE)
    cat(format4Bugs(data), file = dataFile, fill = TRUE)
  }
}


genInitsFile <- function(n.chains, inits, initsFileStem) {
  if(length(inits) == n.chains){   #Modified by Marcos
     for (i in 1:n.chains) {
       file <- paste(initsFileStem, i, ".txt", sep="")
       if (is.function(inits)) cat(format4Bugs(inits()), file=file, fill = TRUE)
       else cat(format4Bugs(inits[[i]]), file=file, fill = TRUE) #Modified by Marcos
     }
  }
  else {                           #Modified by Marcos
     for (i in 1:n.chains) {
       file <- paste(initsFileStem, i, ".txt", sep="")
       if (is.function(inits)) cat(format4Bugs(inits()), file=file, fill = TRUE)
       else cat(format4Bugs(inits[[1]]), file=file, fill = TRUE) #Modified by Marcos
     }
  }
}


#### run bugs
runBugs <- function(bugs=system("which OpenBUGS", TRUE),
                    script,
                    n.chains,
                    workingDir,
                    #useWine=FALSE,
                    #wine = Sys.getenv("WINE"),
                    OpenBugs=TRUE,
                    Windows=TRUE, ## Modified by Marcos
                    verbose = TRUE) {
#  BUGS <- Sys.getenv("BUGS")
#  if (!file.exists(BUGS)) stop(paste(BUGS, "does not exists."))
 # if (!OpenBugs) {
 # if (Windows || useWine) { ## Modified by Marcos
  if(length(n.chains) == 0) stop("The number of MCMC chains in the script file must be specified")
  if(length(workingDir) == 0) stop("A working directory must be specified")
  if (Windows) { ## Modified by Marcos
    if (is.na(pmatch("\"", bugs))) bugs <- paste("\"", bugs, "\"", sep="")
    if (is.na(pmatch("\"", script))) script <- paste("\"", script, "\"", sep="")
    command <- paste(bugs, "/par", script)
  }
  else {
    log.file <- file.path(workingDir, "log.txt") ##Modified by Marcos
    #file.create(log.file,showWarnings = FALSE) ##Modified by Marcos
    command <- paste(bugs, "<", script, ">", log.file)
  }
#  if (useWine) {
#    command <- paste(wine, command)
#
#    ## put a "q" to quit from wine debugg
#    q.tmp <- tempfile("q")
#    on.exit(unlink(q.tmp))
#    cat("q\n", file=q.tmp)
#    command <- paste(command, "< ", q.tmp)
#
#    ## redirect the erorr/warning message of Wine
#    wine.warn <- tempfile("warn")
#    on.exit(unlink(wine.warn))
#    command <- paste(command, ">", wine.warn, " 2>&1 ")
#  }
   
  ## clean up previous coda files
  fnames <- getCodaFileNames(n.chains, workingDir, OpenBugs)
  coda.files <- c(fnames$codaIndexFile, fnames$codaFiles)
##   coda.files <- paste ("coda", 1:n.chains, ".txt", sep="")
##   coda.files <- c("codaIndex.txt", coda.files)
##   coda.files <- file.path(workingDir, coda.files)
  for (i in coda.files) {
    ## cat ("Bugs did not run correctly.\n", file=coda.files[i], append=FALSE)
   if (file.exists(i)) file.remove(i)
  }
  log.file <- file.path(workingDir, "log.txt") ##Modified by Marcos
  if (file.exists(log.file)) file.remove(log.file)
  
  ## execute it!
  cont <- 0
  err <- -1
  while((cont < 2) & (err != 0)){
    cont <- cont + 1
    err <- system(command)
  }
  #if (err == -1) stop("System call to BUGS failed.")
  if (err != 0) stop("System call to BUGS failed.") ## Modified by Marcos
  ## show log
  if (verbose) file.show(log.file) ##Modified by Marcos

  if (!file.exists(coda.files[1])){
    if (Windows) warning("Number of iterations may be to small") ## Modified by Marcos
    stop("BUGS stopped before getting to coda.")
  }
}


#### functions to get the output
getCodaFileNames <- function(n.chains, workingDir, OpenBugs) {
  CODA <- if (OpenBugs) "codaCODA" else "coda"
  INDEX <- if (OpenBugs) "index" else "Index"
  CHAIN <- if (OpenBugs) "chain" else NULL
  coda  <- file.path(workingDir, CODA)
  codaFiles <- paste(coda, CHAIN, 1:n.chains, ".txt", sep="")
  codaIndexFile <- paste(coda, INDEX, ".txt", sep="")
  list(codaFiles=codaFiles, codaIndexFile=codaIndexFile)
}


getBugsOutput <- function(n.chains, workingDir, OpenBugs=TRUE) {
  fnames <- getCodaFileNames(n.chains, workingDir, OpenBugs)
  codaFiles <- fnames$codaFiles
  codaIndexFile <- fnames$codaIndexFile
  if (OpenBugs)  sep <- " "  else sep <- "\t"
  codaIndex <- read.table(codaIndexFile, header=FALSE,
                          sep=sep, as.is=TRUE)

  n.keep <- codaIndex[1, 3] - codaIndex[1, 2] + 1
  ## If data is saved using the OpenBugs interface it may come as "\t" sep
  if(length(n.keep) == 0){                                 ##Modified by Marcos
    sep <- "\t"                                            ##Modified by Marcos
    codaIndex <- read.table(codaIndexFile, header=FALSE,   ##Modified by Marcos
                          sep=sep, as.is=TRUE)             ##Modified by Marcos
    n.keep <- codaIndex[1, 3] - codaIndex[1, 2] + 1        ##Modified by Marcos
  }
  nodes <- codaIndex[, 1]
  n.param <- length(nodes)
  output <- list()
  chain.name <- "chain1"         ##Modified by Marcos
  if (n.chains > 1){             ##Modified by Marcos
   for (i in 2:n.chains) {       ##Modified by Marcos
   k <- paste("chain",i,sep="")  ##Modified by Marcos
   chain.name <- c(chain.name,k) ##Modified by Marcos
   }                             ##Modified by Marcos
  }                              ##Modified by Marcos
  for (i in 1:n.chains) {
    foo <- read.table(codaFiles[i], header=FALSE)
#    iter <- foo[1:n.keep, 1]
    vals <- matrix(foo[,2], n.keep, n.param)
#    dat <- as.data.frame(cbind(iter, vals))
     dat <- as.data.frame(vals)
#    names(dat) <- c("iter", nodes)
    names(dat) <- nodes
    output[[chain.name[i]]] <- dat
  }    
  
  output[[i+1]] <- n.keep  ##Modified by Marcos
  output[[i+2]] <- n.keep*n.chains ##Modified by Marcos

  all <- list (n.chains=n.chains, n.iter=NA, n.burnin=NA, ##Modified by Marcos
               n.thin=NA, n.keep=output[[i+1]], n.sims=output[[i+2]]) ##Modified by Marcos
  output <- output[-(n.chains+2)] ##Modified by Marcos
  output <- output[-(n.chains+1)] ##Modified by Marcos
  all <- c(all, output) ##Modified by Marcos 
  class(all) <- "rbugs" ##Modified by Marcos

  all
  #output
}

jags2 <- function (data, inits, parameters.to.save, model.file = "model.bug",
  n.chains = 3, n.iter = 2000, n.burnin = floor(n.iter/2),
  n.thin = max(1, floor((n.iter - n.burnin)/1000)),
  DIC = TRUE, jags.path = "", working.directory = NULL, clearWD = TRUE,
  refresh = n.iter/50)
{
  if (!is.null(working.directory)) {
    working.directory <- path.expand(working.directory)
    savedWD <- getwd()
    setwd(working.directory)
    on.exit(setwd(savedWD))
  }
  else {
    savedWD <- getwd()
    working.directory <- savedWD
  }

  redo <- ceiling(n.iter - n.burnin)

 # data.list <- lapply(as.list(data), get, pos = parent.frame(2))
#  names(data.list) <- as.list(data)


  if(is.list(data)){
    data.list <- data
    lapply(data.list, dump, append=TRUE, file="jagsdata.txt", envir=parent.frame(1))
  }
  else{
    if (!(length(data) == 1 && is.vector(data) && is.character(data) &&
          (regexpr("\\.txt$", data) > 0))) {
      data.list <- lapply(as.list(data), get, pos = parent.frame(1))
      names(data.list) <- as.list(data)
    }
    else {
      if (all(basename(data) == data)) {
        try(file.copy(file.path(working.directory, data), data, overwrite = TRUE))
      }
      if (!file.exists(data)) {
        stop("File", data, "does not exist.")
      }
      data.list <- data
    }

  }
  lapply(names(data.list), dump, append=TRUE, file="jagsdata.txt")
  data <- read.jagsdata("jagsdata.txt")

  if (is.function(model.file)) {
    temp <- tempfile("model")
    temp <- if (is.R() || .Platform$OS.type != "windows") {
        paste(temp, "txt", sep = ".")
    }
    else {
        gsub("\\.tmp$", ".txt", temp)
    }
    write.model(model.file, con = temp)
    model.file <- gsub("\\\\", "/", temp)
    if (!is.R())
        on.exit(file.remove(model.file), add = TRUE)
  }
  jags.call <- if (jags.path == "") {
    "jags"
  }
  else {
    jags.path <- win2unixdir(jags.path)
    paste(jags.path, "jags", sep = "")
  }
  no.inits <- FALSE
  inits.files <- NULL
  chain.names <- NULL
  if(is.null(inits)){
    no.inits <- TRUE
  }
  else if (is.function(inits)){
    for (i in 1:n.chains) {
      initial.values <- inits()
      inits.files <- c(inits.files, paste("jagsinits", i, ".txt", sep = ""))
      chain.names <- c(chain.names, paste("chain(", i, ")", sep = ""))
      with(initial.values, dump(names(initial.values), file = paste("jagsinits", i, ".txt", sep = "")))
    }
  }
  else if (is.list(inits)){
    if (length(inits)==n.chains){
      for (i in 1:n.chains) {
        initial.values <- inits[[i]]
        inits.files <- c(inits.files, paste("jagsinits", i, ".txt", sep = ""))
        chain.name <- c(chain.names, paste("chain(", i, ")", sep = ""))
        with(initial.values, dump(names(initial.values), file = paste("jagsinits", i, ".txt", sep = "")))
      }
    } else {
      stop(message="initial value must be specified for all of chains")
    }
  }

  if (DIC){
    parameters.to.save <- c(parameters.to.save, "deviance")
    #load.module("dic", quiet = TRUE)
  }


  cat("model clear\ndata clear\n",
      if(DIC){
        "load dic\n"
      },
      "model in ", "\"", model.file, "\"", "\n",
      "cd ", "\"", working.directory, "\"", "\n",
      "data in ", "\"jagsdata.txt\"", "\n",
      "compile, nchains(", n.chains, ")", "\n",
      if(!no.inits){
        paste("inits in \"", inits.files, "\", ", chain.names, "\n", sep = "")
      },
      "initialize", "\n",
      "update ", n.burnin, ", by(", refresh, ")\n",
      paste("monitor ", parameters.to.save, ", thin(", n.thin, ")\n", sep = ""),
      "update ", redo, ", by(", refresh, ")\n",
      "coda *\n", sep = "", file = "jagsscript.txt")

  system(paste(jags.call, "jagsscript.txt"))

  fit <- jags.sims(parameters.to.save = parameters.to.save,
      n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin,
      n.thin = n.thin, DIC = DIC)

  if (clearWD) {
      file.remove(c("jagsdata.txt", "CODAindex.txt", inits.files,
          "jagsscript.txt", paste("CODAchain", 1:n.chains,
              ".txt", sep = "")))
  }
  class(fit) <- "bugs"
  return(fit)
}

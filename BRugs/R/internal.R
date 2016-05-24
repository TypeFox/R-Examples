### Functions to run a single OpenBUGS API command string

.BugsCmd <- function(command)
{
    unlist(.OpenBUGS(command, "BugsCmd"))
}

.CmdInterpreter <- function(command)
{
    unlist(.OpenBUGS(command, "CmdInterpreter"))
}

.Integer <- function(command)
{
    unlist(.OpenBUGS(command, "Integer"))
}

.CharArray <-  function(command, arg)
{
    unlist(.OpenBUGS(command, "CharArray", arg))
}

.RealArray <-  function(command, arg)
{
    unlist(.OpenBUGS(command, "RealArray", arg))
}


.OpenBUGS.cmdtypes <- c("CmdInterpreter","Integer","CharArray","RealArray","BugsCmd")

.OpenBUGS <- function(cmds, cmdtypes=NULL, args=NULL) {
    ncmds <- length(cmds)
    if (is.null(cmdtypes)) cmdtypes <- rep("CmdInterpreter", ncmds)
    if (is.null(args)) args <- as.list(rep(NA, ncmds))
    stopifnot(ncmds==length(cmdtypes))
    stopifnot(ncmds==length(args))
    .OpenBUGS.platform(cmds, cmdtypes, args)
}

dquote <- function(x){
    paste("\"", x, "\"", sep="")
}

.OpenBUGS.helper <- function(cmds, cmdtypes, args) {
    ncmds <- length(cmds)
    if (ncmds > 99999) stop("Maximum number of OpenBUGS API commands exceeded")
    tempDir <- getOption("BRugsTmpdir")
    ## Don't want internalize/externalize to overwrite the command
    ## output buffer, so redirect its output to a separate trash can.
    trashDir <- file.path(tempDir, "trash", fsep="/")
    extFile <- getOption("BRugsExtFile")
    cmdFile <- paste(tempDir, "cmds.txt", sep="/")
    bugsPath <- system.file("exec", paste("BugsHelper", if(.Platform$OS.type == "windows") ".exe", sep=""), package="BRugs")
    shcmd <- paste(dquote(bugsPath), dquote(tempDir), dquote(trashDir), dquote(extFile), dquote(cmdFile), dquote(ncmds))
    for (i in 1:ncmds) {
        if (cmdtypes[i] %in% c("CharArray","RealArray"))
            cat(args[[i]], file=paste(tempDir, "/input",i,".txt", sep=""))
    }
    cmd.id <- match(cmdtypes, .OpenBUGS.cmdtypes) - 1
    write(rbind(cmds, cmd.id), cmdFile)
    res <- system(shcmd)
    handleRes(res)
    out <- vector(ncmds, mode="list")
    for (i in seq_along(cmds)){
        if (cmdtypes[i] %in% c("Integer","CharArray","RealArray"))
            out[[i]] <- scan(paste(tempDir,"/output",i,".txt",sep=""),
                             switch(cmdtypes[i],
                                    "Integer" = integer(),
                                    "CharArray" = character(),
                                    "RealArray" = double()),
                             quiet=TRUE)
    }
    out
}

handleRes <- function(res)
{
    maintainer <- maintainer("BRugs")
    errors <- c("Internal \"trap\" error in OpenBUGS, or non-existent module or procedure called.",
                "An OpenBUGS procedure was called with the wrong type of argument.",
                "An OpenBUGS procedure was called with the wrong signature.")
    ## If a library call ends in a trap, then error code 1 will be returned from BugsHelper on Linux
    ## On Windows it shouldn't even get this far after a trap. TODO see if the trap message is shown.
    if (res > 0) {
        buf <- readLines(file.path(tempdir(), "buffer.txt"))
        trap <- grep("Sorry something went wrong", buf, value=TRUE)
        if(length(trap) > 0) message(trap[1])
        stop(errors[res])#, "\nPlease report this bug to ", maintainer)
    }
}

.SamplesGlobalsCmd <- function(node){
    options.old <- options()
    options(scipen=20) # don't pass numbers in scientific notation to OpenBUGS
    commands <- c(paste("SamplesEmbed.beg :=", getOption("BRugsSamplesBeg")),
                  paste("SamplesEmbed.end :=", getOption("BRugsSamplesEnd")),
                  paste("SamplesEmbed.firstChain :=", getOption("BRugsSamplesFirstChain")),
                  paste("SamplesEmbed.lastChain :=", getOption("BRugsSamplesLastChain")),
                  paste("SamplesEmbed.thin :=", getOption("BRugsSamplesThin")),
                  paste("SamplesEmbed.SetVariable(", sQuote(node), ")", sep=""),
                  paste("BugsMappers.SetPrec(", getOption("BRugsPrec"), ")", sep="")
                  )
    options(options.old)
    paste(commands, collapse="; ")
}

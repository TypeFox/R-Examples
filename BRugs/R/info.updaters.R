
## display updaters sorted by node name
infoUpdatersbyName <- function()
{
    command <- "BugsEmbed.NotCompiledGuard; BugsEmbed.UpdatersByName"
    .CmdInterpreter(command)
    buffer <- file.path(tempdir(), "buffer.txt")
    if (readLines(buffer)[1]=="BugsCmds:NotCompiled")
        stop("Model not compiled")
    buffer <- file.path(tempdir(), "Updater types.txt")
    result <- read.table(buffer, sep="\t", skip=1, as.is=TRUE,
                       row.names=2, col.names=c("empty", "Node", "Type", "Size", "Depth"))[,-1]
    ## strip leading and trailing spaces
    for (i in 1:2) { 
        result[,i] <- gsub("^ +", "\\1", result[,i])
        result[,i] <- gsub(" +$", "\\1", result[,i])
    } 
    rownames(result) <- gsub("^ +", "", rownames(result))
    rownames(result) <- gsub(" +$", "", rownames(result))
    unlink(buffer)
    result
}

## display updaters sorted by node depth in graph
infoUpdatersbyDepth <- function()
{
    command <- "BugsEmbed.NotCompiledGuard; BugsEmbed.UpdatersByDepth"
    .CmdInterpreter(command)
    buffer <- file.path(tempdir(), "buffer.txt")
    if (readLines(buffer)[1]=="BugsCmds:NotCompiled")
        stop("Model not compiled")
    buffer <- file.path(tempdir(), "Updater types.txt")
    result <- read.table(buffer, sep="\t", skip=1, as.is=TRUE,
                       row.names=2, col.names=c("empty", "Node", "Type", "Size", "Depth"))[,-1]
    ## strip leading and trailing spaces
    for (i in 1:2) { 
        result[,i] <- gsub("^ +", "\\1", result[,i])
        result[,i] <- gsub(" +$", "\\1", result[,i])
    } 
    rownames(result) <- gsub("^ +", "", rownames(result))
    rownames(result) <- gsub(" +$", "", rownames(result))
    unlink(buffer)
    result
}

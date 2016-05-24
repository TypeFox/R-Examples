"modelInits" <-
function(fileName, chainNum = NULL)
#   Load initial values for OpenBUGS model
{
    if(is.null(chainNum))
        chainNum <- getChain() + seq(along = fileName) - 1
    if(!is.numeric(chainNum))
        stop("chainNum ", "must be numeric")
    if(length(fileName) != length(chainNum))
        stop("length(chainNum) ", "must be equal to the number of filenames given")
    chainNum <- as.integer(chainNum)
    path <- dirname(fileName)
    path <- ifelse(path == ".", getwd(), path)
    fileName <- file.path(path, basename(fileName))
    fileExist <- !file.exists(fileName)
    if(any(fileExist))
        stop("File(s) ", fileName[fileExist], " do(es) not exist.")
    for(i in seq(along = fileName)){        
        if(file.info(fileName[i])$isdir) 
            stop(fileName[i], " is a directory, but a file is required.")
        filename <- sQuote(fileName[i])
        command <- paste("BugsEmbed.SetFilePath(", filename, 
        "); BugsEmbed.LoadInitsGuard; BugsEmbed.chain := ",
        as.character(chainNum[i]), "; BugsEmbed.LoadInits")
        if (!is.R()){
          command <- gsub ("\\\\", "/", command)
            command <- gsub ("//", "/", command)
        }
        .CmdInterpreter(command)
        if(getOption("BRugsVerbose")){
            message("Initializing chain ", chainNum[i], ": ", sep="")
            buffer()
        }
        options("BRugsNextChain" = chainNum[i] + 1)        
    }
    invisible()
}

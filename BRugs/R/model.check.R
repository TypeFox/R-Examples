"modelCheck" <-
function(fileName)
#   Check that OpenBUGS model is syntactically correct
{
    path <- dirname(fileName)
    path <- if(path == ".") getwd() else path
    fileName <- file.path(path, basename(fileName))
    if(!file.exists(fileName))
        stop("File ", fileName, " does not exist")
    if(file.info(fileName)$isdir) 
        stop(fileName, " is a directory, but a file is required")
    command <- paste("BugsEmbed.SetFilePath(", sQuote(fileName), 
        ");BugsEmbed.ParseGuard;BugsEmbed.Parse", sep = "")
    if (!is.R()) {
        command <- gsub ("\\\\", "/", command)
        command <- gsub ("//", "/", command)
    }
    .CmdInterpreter(command)
    .initGlobals()
    if(getOption("BRugsVerbose")) 
        buffer()
}

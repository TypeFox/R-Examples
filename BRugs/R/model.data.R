"modelData" <-
function(fileName = "data.txt")
{
#   Load data for OpenBUGS model
    for(i in fileName){
        path <- dirname(i)
        path <- if(path == ".") getwd() else path
        fileNm <- file.path(path, basename(i))
        if(!file.exists(fileNm))
            stop("File ", fileNm, " does not exist")
        if(file.info(fileNm)$isdir) 
            stop(fileNm, " is a directory, but a file is required")
        command <- paste("BugsEmbed.SetFilePath(", sQuote(fileNm), 
            ");BugsEmbed.LoadDataGuard;BugsEmbed.LoadData", sep = "")
        if (!is.R()){
          command <- gsub ("\\\\", "/", command)
            command <- gsub ("//", "/", command)
        }
        .CmdInterpreter(command)
        if(getOption("BRugsVerbose"))     
            buffer()
    }
}

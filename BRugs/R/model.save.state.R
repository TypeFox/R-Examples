"modelSaveState" <- function(stem)
{
## Saves the sate of each chain in OpenBUGS model
    if(!is.character(stem) || length(stem)!=1)
        stop("'stem' must be character of length 1")
    if(dirname(stem) == ".") 
        stem <- file.path(getwd(), basename(stem))
    command <- paste("BugsEmbed.UpdateGuard", ";BugsEmbed.WriteChains(", sQuote(stem), ")")
    .CmdInterpreter(command)
    if(getOption("BRugsVerbose")) 
      buffer()
}

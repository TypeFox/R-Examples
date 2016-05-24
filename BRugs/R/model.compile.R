"modelCompile" <-
function(numChains = 1)
#   Compile OpenBUGS model
{
    if(!is.numeric(numChains))
        stop("numChains ", "must be numeric")
    numChains <- as.integer(numChains)
    command <- paste("BugsEmbed.CompileGuard",              
        ";BugsEmbed.numChains :=", as.character(numChains), "; BugsEmbed.Compile", sep = "")
    .CmdInterpreter(command)
    samplesSetFirstChain(1)
    samplesSetLastChain(numChains)
    options("BRugsNextChain" = 1)
    if(getOption("BRugsVerbose")) 
        buffer()
}

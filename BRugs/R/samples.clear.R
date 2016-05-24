"samplesClear" <-
function(node)
#   Clear a sample monitor
{
    nodeName <- sQuote(node)
    command <- paste("SamplesEmbed.SetVariable(", nodeName, 
                            ");SamplesEmbed.HistoryGuard;SamplesEmbed.Clear")
    .CmdInterpreter(command)
    if(getOption("BRugsVerbose")) 
        buffer()
}

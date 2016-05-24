"ranksSet" <-
function(node)
#   Set a ranks monitor for vector quantity node in OpenBUGS model
{
    nodeName <- sQuote(node)
    command <- paste("RanksEmbed.SetVariable(", nodeName, "); RanksEmbed.SetGuard;",
                             "RanksEmbed.Set")
    .CmdInterpreter(command)
    if(getOption("BRugsVerbose")) 
        buffer()
}

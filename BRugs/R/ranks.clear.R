"ranksClear" <-
function(node)
#   Clears a ranks monitor for vector quantity in OpenBUGS model
{
    nodeName <- sQuote(node)
    command <- paste("RanksEmbed.SetVariable(", nodeName, "); RanksEmbed.StatsGuard;",
                                 "RanksEmbed.Clear")
    .CmdInterpreter(command)
    if(getOption("BRugsVerbose")) 
        buffer()
}

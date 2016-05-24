"summaryClear" <-
function(node)
#   Clear summary monitor for node in WinBUGS model
{
    nodeName <- sQuote(node)
    for(i in seq(along=nodeName)){
        command <- paste("SummaryEmbed.SetVariable(", nodeName[i], "); SummaryEmbed.StatsGuard;",
                         "SummaryEmbed.Clear")
        .CmdInterpreter(command)
        buffer <- file.path(tempdir(), "buffer.txt")
        rlb <- readLines(buffer)
        if(getOption("BRugsVerbose"))
            message("Variable ", nodeName[i], ": ", rlb)
    }
    invisible()
}
